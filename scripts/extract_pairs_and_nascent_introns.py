import os

os.environ[
    'OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

import argparse
import json
import sys

import pysam
from interval import interval as py_interval
import pandas as pd
from typing import NamedTuple, Optional

import numpy as np
from pathlib import Path
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.DEBUG,
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout)])


class GenomicRange(NamedTuple):
    chromosome: str
    start: int
    end: int
    strand: Optional[str] = None

    def unstranded_bed_string(self):
        return f"{self.chromosome}\t{self.start}\t{self.end}\n"


class ChromsomeAndStrand(NamedTuple):
    chromosome: str
    strand: str


class IntronsIndex:

    def __init__(self, introns_bed_file: Path, fai_index_file: Path) -> None:
        fai_df = pd.read_csv(fai_index_file, sep='\t',
                             names=['chromosome', 'length', 'offset', 'linebases', 'linewidth'])
        fai_df['chromosome'] = fai_df['chromosome'].astype(str)
        self.index_by_chrom_and_strand: dict[ChromsomeAndStrand, dict] = {}
        for strand in ('+', '-'):
            for chromosome in fai_df['chromosome']:
                self.index_by_chrom_and_strand[ChromsomeAndStrand(chromosome=chromosome, strand=strand)] = {}

        introns_df = pd.read_csv(introns_bed_file, sep='\t',
                                 names=['chromosome', 'start', 'end', 'name', 'score', 'strand'])
        introns_df['chromosome'] = introns_df['chromosome'].astype(str)
        logging.info(f"Loading introns for indexing")
        for _, row in introns_df.iterrows():
            intron = GenomicRange(chromosome=row['chromosome'], start=row['start'], end=row['end'],
                                  strand=row['strand'])
            chrom_and_strand = ChromsomeAndStrand(chromosome=row['chromosome'], strand=row['strand'])
            self.index_by_chrom_and_strand[chrom_and_strand][
                len(self.index_by_chrom_and_strand[chrom_and_strand]) + 1] = intron
        largest_index = max([max(numeric_index.keys()) if numeric_index.keys() else 0
                             for numeric_index in self.index_by_chrom_and_strand.values()])
        selected_int_type = np.uint16 if largest_index <= np.iinfo(np.uint16).max else np.uint32

        self.genomic_index: dict[ChromsomeAndStrand, np.array] = {}
        for strand in ('+', '-'):
            for chromosome, length in zip(fai_df['chromosome'], fai_df['length']):
                self.genomic_index[ChromsomeAndStrand(chromosome=chromosome, strand=strand)] = np.zeros(length,
                                                                                                        dtype=selected_int_type)
        for chromosome_and_strand, numeric_index in self.index_by_chrom_and_strand.items():
            for index, intron in numeric_index.items():
                self.genomic_index[chromosome_and_strand][intron.start:intron.end] = index

    def find_overlapping_intron(self, chromosome: str, strand: str, position: int) -> Optional[GenomicRange]:
        chrom_and_strand = ChromsomeAndStrand(chromosome=chromosome, strand=strand)
        index = self.genomic_index[chrom_and_strand][position]
        return self.index_by_chrom_and_strand[chrom_and_strand][index] if index != 0 else None


def read_is_in_forward_pair(read: pysam.AlignedSegment, strandendess_type: str) -> bool:
    assert strandendess_type in ['1', '2'], "strandendess_type must be either '1' or '2'"
    if strandendess_type == '1':
        if (read.is_read1 and read.is_forward) or (read.is_read2 and read.is_reverse):
            return True
        else:
            return False
    elif strandendess_type == '2':
        if (read.is_read2 and read.is_forward) or (read.is_read1 and read.is_reverse):
            return True
        else:
            return False


def extract_id_of_invalid_reads(bamfile_input_path: Path) -> set[str]:
    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")
    invalid_ids: set[str] = set()
    for i, read in enumerate(bamfile_input):
        if i % 1_000_000 == 0:
            logging.info(f"Finding invalid reads: {i} reads")
        if read.is_secondary:
            invalid_ids.add(read.query_name)
    bamfile_input.close()
    return invalid_ids


def extract_and_save_unique_pairs(input_folder: Path,
                                  output_folder: Path,
                                  strandendess_type: str,
                                  introns_bed_file: Path,
                                  fai_index_file: Path,
                                  bam_file_name="Aligned.sortedByCoord.out.bam"
                                  ) -> None:
    assert strandendess_type in ['1', '2'], "strandendess_type must be either '1' or '2'"

    bamfile_input_path = input_folder / bam_file_name

    output_bam_file_forward_path = output_folder / 'forward.bam'
    output_bam_file_reverse_path = output_folder / 'reverse.bam'

    output_bed_forward_pairs = output_folder / 'forward_pairs.bed'
    output_bed_reverse_pairs = output_folder / 'reverse_pairs.bed'

    output_bed_forward_nascent_introns = output_folder / 'forward_nascent_introns.bed'
    output_bed_reverse_nascent_introns = output_folder / 'reverse_nascent_introns.bed'

    for file_name in [output_bed_forward_pairs, output_bed_reverse_pairs,
                      output_bed_forward_nascent_introns, output_bed_reverse_nascent_introns]:
        open(file_name, 'w').close()  # Create empty files to append on

    output_json_read_count_file = output_folder / 'read_counts.json'

    invalid_ids = extract_id_of_invalid_reads(bamfile_input_path)

    reads_1: dict[str, pysam.AlignedSegment] = {}
    reads_2: dict[str, pysam.AlignedSegment] = {}

    intervals_forward_pairs: list[GenomicRange] = []
    intervals_reverse_pairs: list[GenomicRange] = []

    intervals_forward_nascent_introns: list[GenomicRange] = []
    intervals_reverse_nascent_introns: list[GenomicRange] = []

    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")

    introns_index = IntronsIndex(introns_bed_file, fai_index_file)

    valid_reads_forward = 0
    valid_reads_reverse = 0
    for i, read in enumerate(bamfile_input):
        if i % 1_000_000 == 0:
            logging.info(f"Computing covered intervals: {i} reads")
            for output_file, intervals in [(output_bed_forward_pairs, intervals_forward_pairs),
                                           (output_bed_reverse_pairs, intervals_reverse_pairs),
                                           (output_bed_forward_nascent_introns, intervals_forward_nascent_introns),
                                           (output_bed_reverse_nascent_introns, intervals_reverse_nascent_introns)]:
                with open(output_file, 'a') as file:
                    while intervals:
                        file.write(intervals.pop().unstranded_bed_string())

        if read.query_name in invalid_ids:
            continue

        if read.is_read1 and read.query_name not in reads_2:
            reads_1[read.query_name] = read
            continue
        elif read.is_read2 and read.query_name not in reads_1:
            reads_2[read.query_name] = read
            continue

        if read.is_read1 and read.query_name in reads_2:
            read_1 = read
            read_2 = reads_2.pop(read.query_name)
        elif read.is_read2 and read.query_name in reads_1:
            read_1 = reads_1.pop(read.query_name)
            read_2 = read
        else:
            assert False

        if read_1.reference_name != read_2.reference_name:
            invalid_ids.add(read_1.query_name)
            continue

        interval_union = py_interval(*(read_1.get_blocks() + read_2.get_blocks()))
        interval_union = sorted(list(interval_union), key=lambda x: x[0])

        if read_is_in_forward_pair(read=read_1, strandendess_type=strandendess_type):
            intervals_forward_pairs.extend([GenomicRange(chromosome=read_1.reference_name,
                                                         start=int(x[0]),
                                                         end=int(x[1]),
                                                         strand='+') for x in interval_union])
            valid_reads_forward += 2
            # suspected_polymerase_position is inclusive (we expect pol-II to be located there).
            # Intervals in BED format have left part (start) inclusive and right part (end) exclusive, which
            # is why we are adjusting by 1.
            suspected_polymerase_position = max([int(x[1]) for x in interval_union]) - 1
            overlapping_intron = introns_index.find_overlapping_intron(chromosome=read_1.reference_name,
                                                                       strand='+',
                                                                       position=suspected_polymerase_position)
            if overlapping_intron is not None:
                intervals_forward_nascent_introns.append(GenomicRange(chromosome=read_1.reference_name,
                                                                      start=overlapping_intron.start,
                                                                      end=suspected_polymerase_position + 1,
                                                                      strand='+'))

        else:
            intervals_reverse_pairs.extend([GenomicRange(chromosome=read_1.reference_name,
                                                         start=int(x[0]),
                                                         end=int(x[1]),
                                                         strand='-') for x in interval_union])
            valid_reads_reverse += 2

            suspected_polymerase_position = min([int(x[0]) for x in interval_union])
            overlapping_intron = introns_index.find_overlapping_intron(chromosome=read_1.reference_name,
                                                                       strand='-',
                                                                       position=suspected_polymerase_position)
            if overlapping_intron is not None:
                intervals_reverse_nascent_introns.append(GenomicRange(chromosome=read_1.reference_name,
                                                                      start=suspected_polymerase_position,
                                                                      end=overlapping_intron.end,
                                                                      strand='-'))
    for output_file, intervals in [(output_bed_forward_pairs, intervals_forward_pairs),
                                   (output_bed_reverse_pairs, intervals_reverse_pairs),
                                   (output_bed_forward_nascent_introns, intervals_forward_nascent_introns),
                                   (output_bed_reverse_nascent_introns, intervals_reverse_nascent_introns)]:
        with open(output_file, 'a') as file:
            while intervals:
                file.write(intervals.pop().unstranded_bed_string())
    bamfile_input.close()

    with open(output_json_read_count_file, 'w') as output_json_file:
        json.dump({'selected_reads_forward': valid_reads_forward,
                   'selected_reads_reverse': valid_reads_reverse},
                  output_json_file)

    # Write .bam files
    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")
    bamfile_output_forward = pysam.AlignmentFile(output_bam_file_forward_path, "wb",
                                                 template=bamfile_input)
    bamfile_output_reverse = pysam.AlignmentFile(output_bam_file_reverse_path, "wb",
                                                 template=bamfile_input)
    for i, read in enumerate(bamfile_input):
        if i % 1_000_000 == 0:
            logging.info(f"Writing reads: {i} reads")
        if read.query_name in invalid_ids:
            continue

        if read_is_in_forward_pair(read=read, strandendess_type=strandendess_type):
            bamfile_output_forward.write(read)
        else:
            bamfile_output_reverse.write(read)

    bamfile_output_forward.close()
    bamfile_output_reverse.close()
    bamfile_input.close()

    logging.info("Sorting and indexing output .bam files")
    for bam_file_path in (output_bam_file_forward_path, output_bam_file_reverse_path):
        tmp_file_path = bam_file_path.parent / 'tmp_sorted.bam'
        pysam.sort("-o", str(tmp_file_path), str(bam_file_path))
        os.rename(tmp_file_path, bam_file_path)
        pysam.index(str(bam_file_path))
    logging.info("Extracting pairs and nascent introns finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder')
    parser.add_argument('--output_folder')
    parser.add_argument('--strandendess_type')
    parser.add_argument('--introns_bed_file')
    parser.add_argument('--fai_index_file')
    args = parser.parse_args()
    extract_and_save_unique_pairs(input_folder=Path(args.input_folder),
                                  output_folder=Path(args.output_folder),
                                  strandendess_type=args.strandendess_type,
                                  introns_bed_file=Path(args.introns_bed_file),
                                  fai_index_file=Path(args.fai_index_file))
