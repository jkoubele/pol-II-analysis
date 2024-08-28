import argparse
from pathlib import Path

import pandas as pd
from pybedtools import BedTool, Interval
from tqdm import tqdm


def extract_genomic_features(genome_folder: Path, gtf_file_name: str, gtf_source: str) -> None:
    assert gtf_source in ('ensembl', 'gencode'), "gtf_source must be either 'ensembl' or 'gencode'."
    gtf_df = pd.read_csv(genome_folder / gtf_file_name,
                         header=4,
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                                'attribute'],
                         delimiter="\t")

    print(f"Possible features: {gtf_df['feature'].unique()}")

    # Processing of the reference genome to obtain introns
    gtf_records = BedTool(genome_folder / gtf_file_name)
    genes = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting genes') if record.fields[2] == 'gene']).sort()
    exons = BedTool(
        [record for record in tqdm(gtf_records, desc='Extracting exons') if record.fields[2] == 'exon']).sort()

    if gtf_source == 'ensembl':
        utr_3_prime = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting three_prime_utr") if
             record.fields[2] == 'three_prime_utr']).sort()
        utr_5_prime = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting five_prime_utr") if
             record.fields[2] == 'five_prime_utr']).sort()
        utr_3_prime.saveas(genome_folder / 'utr_3_prime.bed')
        utr_5_prime.saveas(genome_folder / 'utr_5_prime.bed')
        introns = genes.subtract(exons, s=True).subtract(utr_3_prime, s=True).subtract(utr_5_prime, s=True).sort()
    elif gtf_source == 'gencode':
        utr = BedTool(
            [record for record in tqdm(gtf_records, desc="Extracting UTR") if
             record.fields[2] == 'UTR']).sort()
        utr.saveas(genome_folder / 'utr.bed')
        introns = genes.subtract(exons, s=True).subtract(utr, s=True).sort()
    else:
        assert False, "Invalid gtf_source."

    introns = BedTool([Interval(chrom=x.chrom, start=x.start, end=x.end,
                                name=x.fields[-1].split()[1][1:-2], score='.',
                                strand=x.strand)
                       for x in tqdm(introns)])
    introns = introns.merge(s=True, c=[4, 5, 6], o='distinct').sort()

    genes.saveas(genome_folder / 'genes.bed')
    exons.saveas(genome_folder / 'exons.bed')
    introns.saveas(genome_folder / 'introns.bed')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome_folder')
    parser.add_argument('--gtf_file_name')
    parser.add_argument('--gtf_source')
    args = parser.parse_args()
    extract_genomic_features(genome_folder=Path(args.genome_folder),
                             gtf_file_name=args.gtf_file_name,
                             gtf_source=args.gtf_source)
