import argparse
import json
from collections import defaultdict
from pathlib import Path
from typing import Optional

import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--project_folder',
                        default='/cellfile/datapublic/jkoubele/celegans_mutants')
    args = parser.parse_args()
    project_folder = Path(args.project_folder)

    sample_annotation = pd.read_csv(project_folder / 'sample_annotation' / 'sample_annotation.tsv', sep='\t')
    library_sizes: list[int] = []
    for sample_name in sample_annotation['sample_name']:
        with open(project_folder / 'coverage' / sample_name / 'read_counts.json') as file:
            read_counts_json = json.load(file)
        library_sizes.append(read_counts_json['selected_reads_forward'] + read_counts_json['selected_reads_reverse'])

    sample_annotation['library_size'] = library_sizes
    sample_annotation.to_csv(project_folder / 'sample_annotation' / 'sample_annotation_with_library_size.tsv',
                             sep='\t', index=False)

    sample_annotation = sample_annotation.set_index('sample_name', drop=False)
    intron_read_counts: Optional[pd.DataFrame] = None

    for sample_name in sample_annotation['sample_name']:
        slopes_df = pd.read_csv(project_folder / 'intron_slopes' / sample_name / 'slopes_by_definition.tsv',
                                sep='\t')
        slopes_df = slopes_df.dropna()
        if intron_read_counts is None:
            intron_read_counts = slopes_df[['name', 'chromosome', 'start', 'end', 'strand', 'length']]
        intron_read_counts[sample_name] = (slopes_df['coverage_5_prime'] - slopes_df['coverage_3_prime'])

    intron_read_counts = intron_read_counts.rename(columns={'name': 'gene'})
    # C. Elegans have only  361 out of 74k introns overlapping --> we just drop them:
    intron_read_counts = intron_read_counts[[',' not in gene for gene in intron_read_counts['gene']]]

    intron_counter = defaultdict(int)
    intron_numbers: list[int] = []
    for gene in intron_read_counts['gene']:
        intron_counter[gene] += 1
        intron_numbers.append(intron_counter[gene])

    intron_read_counts['intron_number'] = intron_numbers
    intron_read_counts['intron_name'] = [f"{gene}_{number}" for gene, number in zip(intron_read_counts['gene'],
                                                                                    intron_read_counts[
                                                                                        'intron_number'])]

    intron_read_counts.insert(1, 'intron_number', intron_read_counts.pop('intron_number'))
    intron_read_counts.insert(1, 'intron_name', intron_read_counts.pop('intron_name'))

    intron_read_counts.to_csv(project_folder / 'aggregated_intronic_counts' / 'intron_read_counts.tsv', sep='\t',
                              index=False)
