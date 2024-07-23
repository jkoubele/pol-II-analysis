import os

os.environ[
    'OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

import argparse
import json
from pathlib import Path

import pysam
from pybedtools import BedTool, Interval
from interval import interval as py_interval


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
            print(f"Finding invalid reads: {i}", flush=True)
        if read.is_secondary:
            invalid_ids.add(read.query_name)
    bamfile_input.close()
    return invalid_ids


def extract_and_save_unique_pairs(input_folder: Path,
                                  output_folder: Path,
                                  strandendess_type: str,
                                  bam_file_name="Aligned.sortedByCoord.out.bam") -> None:
    assert strandendess_type in ['1', '2'], "strandendess_type must be either '1' or '2'"

    bamfile_input_path = input_folder / bam_file_name

    output_bam_file_forward_path = output_folder / 'forward.bam'
    output_bam_file_reverse_path = output_folder / 'reverse.bam'

    output_bed_file_forward_path = output_folder / 'forward.bed.gz'
    output_bed_file_reverse_path = output_folder / 'reverse.bed.gz'

    output_json_read_count_file = output_folder / 'read_counts.json'

    invalid_ids = extract_id_of_invalid_reads(bamfile_input_path)

    reads_1: dict[str, pysam.AlignedSegment] = {}
    reads_2: dict[str, pysam.AlignedSegment] = {}

    intervals_forward: list[Interval] = []
    intervals_reverse: list[Interval] = []

    bamfile_input = pysam.AlignmentFile(bamfile_input_path, "rb")

    valid_reads_forward = 0
    valid_reads_reverse = 0
    for i, read in enumerate(bamfile_input):
        if i % 1_000_000 == 0:
            print(f"Computing covered intervals: {i}", flush=True)
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
            intervals_forward.extend([Interval(chrom=read_1.reference_name,
                                               start=int(x[0]),
                                               end=int(x[1]),
                                               strand='+') for x in interval_union])
            valid_reads_forward += 2
        else:
            intervals_reverse.extend([Interval(chrom=read_1.reference_name,
                                               start=int(x[0]),
                                               end=int(x[1]),
                                               strand='-') for x in interval_union])
            valid_reads_reverse += 2
    bamfile_input.close()

    print("Saving intervals to .bed files", flush=True)
    BedTool(intervals_forward).sort().saveas(output_bed_file_forward_path, compressed=True)
    BedTool(intervals_reverse).sort().saveas(output_bed_file_reverse_path, compressed=True)

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
            print(f"Writing reads: {i}", flush=True)
        if read.query_name in invalid_ids:
            continue

        if read_is_in_forward_pair(read=read, strandendess_type=strandendess_type):
            bamfile_output_forward.write(read)
        else:
            bamfile_output_reverse.write(read)

    bamfile_output_forward.close()
    bamfile_output_reverse.close()
    bamfile_input.close()

    print("Sorting and indexing output .bam files", flush=True)
    for bam_file_path in (output_bam_file_forward_path, output_bam_file_reverse_path):
        tmp_file_path = bam_file_path.parent / 'tmp_sorted.bam'
        pysam.sort("-o", str(tmp_file_path), str(bam_file_path))
        os.rename(tmp_file_path, bam_file_path)
        pysam.index(str(bam_file_path))


if __name__ == "__main__":    
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder')
    parser.add_argument('--output_folder')
    parser.add_argument('--strandendess_type')
    args = parser.parse_args()
    extract_and_save_unique_pairs(input_folder=Path(args.input_folder),
                                  output_folder=Path(args.output_folder),
                                  strandendess_type=args.strandendess_type)
