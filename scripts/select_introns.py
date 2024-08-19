import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
from pybedtools import BedTool, Interval

SJ_FILE_COLUMNS = ['chromosome', 'start', 'end', 'strand', 'intron_motif',
                   'annotated', 'reads_unique', 'reads_multimapped', 'max_overhang']


def load_and_preprocess_slopes(slopes_file: Path, min_r_squared=0.01) -> pd.DataFrame:
    slopes_df = pd.read_csv(slopes_file, sep='\t')
    slopes_df['chromosome'] = slopes_df['chromosome'].astype(str)
    slopes_df = slopes_df.dropna()
    slopes_df = slopes_df[slopes_df['slope'] < 0]
    if 'r_squared' in slopes_df.columns:
        slopes_df = slopes_df[slopes_df['r_squared'] > min_r_squared]
    if 'num_polymerases_per_million_reads' in slopes_df.columns:
        slopes_df = slopes_df[slopes_df['num_polymerases_per_million_reads'] > 0]
    slopes_df = slopes_df.reset_index(drop=True)
    return slopes_df


def load_and_preprocess_sj(sj_file: Path) -> pd.DataFrame:
    sj_df = pd.read_csv(sj_file, sep='\t', names=SJ_FILE_COLUMNS)
    sj_df['chromosome'] = sj_df['chromosome'].astype(str)

    # STAR SJ files are 1-indexed and interval end is inclusive, while .bed files are 0-based and interval end is exclusive.
    # We convert SJ files to .bed-style indexing (note that interval end is unchanged as 0-/1- indexing cancels out with exclusiveness / inclusiveness)
    sj_df['start'] -= 1

    # We omit SJ with unspecified strand    
    sj_df = sj_df[sj_df['strand'] != 0]

    sj_df['strand'] = sj_df['strand'].apply(lambda x: '+' if x == 1 else '-')
    sj_df = sj_df.reset_index(drop=True)
    return sj_df


def select_intron_slopes_by_sj_evidence(slopes_file: Path, sj_file: Path) -> pd.DataFrame:
    max_sj_shift = 5
    min_unique_reads_for_intron_evidence = 5
    min_unique_reads_for_nested_splicing = 1
    nested_splicing_padding = 20

    slopes_df = load_and_preprocess_slopes(slopes_file)
    sj_df = load_and_preprocess_sj(sj_file)
    bedtool_slopes = BedTool([Interval(chrom=row['chromosome'],
                                       start=row['start'],
                                       end=row['end'],
                                       name=index,
                                       strand=row['strand']) for index, row in slopes_df.iterrows()]).sort()

    bedtool_sj = BedTool([Interval(chrom=row['chromosome'],
                                   start=row['start'],
                                   end=row['end'],
                                   name=index,
                                   strand=row['strand']) for index, row in sj_df.iterrows()]).sort()

    intersection_bedtool = bedtool_slopes.intersect(bedtool_sj, sorted=True, s=True, wa=True, wb=True)

    intersections_by_query_id = defaultdict(list)
    for intersection_interval in intersection_bedtool:
        intersections_by_query_id[int(intersection_interval.fields[3])].append(int(intersection_interval.fields[9]))

    evidence_of_sj: list[bool] = []
    evidence_of_nested_sj: list[bool] = []

    for index, row in slopes_df.iterrows():
        overlapping_sj = sj_df.loc[intersections_by_query_id[index]]
        if len(overlapping_sj) == 0:
            evidence_of_sj.append(False)
            evidence_of_nested_sj.append(False)
            continue

        supporting_sj = overlapping_sj[(overlapping_sj['start'] >= row['start'] - max_sj_shift) &
                                       (overlapping_sj['start'] <= row['start'] + max_sj_shift) &
                                       (overlapping_sj['end'] >= row['end'] - max_sj_shift) &
                                       (overlapping_sj['end'] <= row['end'] + max_sj_shift)]

        nested_region_start = row['start'] + nested_splicing_padding
        nested_region_end = row['end'] - nested_splicing_padding

        if nested_region_end <= nested_region_start:
            evidence_of_nested_sj.append(False)
            continue

        evidence_of_sj.append(supporting_sj['reads_unique'].sum() >= min_unique_reads_for_intron_evidence)

        nested_sj = overlapping_sj[((overlapping_sj['start'] >= nested_region_start) &
                                    (overlapping_sj['start'] <= nested_region_end)) |
                                   ((overlapping_sj['end'] >= nested_region_start) &
                                    (overlapping_sj['end'] <= nested_region_end))]
        evidence_of_nested_sj.append(nested_sj['reads_unique'].sum() > min_unique_reads_for_nested_splicing)
    slopes_df['evidence_of_sj'] = evidence_of_sj
    slopes_df['evidence_of_nested_sj'] = evidence_of_nested_sj
    return slopes_df[slopes_df['evidence_of_sj'] & ~slopes_df['evidence_of_nested_sj']]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder_slopes')
    parser.add_argument('--input_folder_sj')
    parser.add_argument('--output_folder')
    args = parser.parse_args()

    selected_slopes_by_definition = select_intron_slopes_by_sj_evidence(
        slopes_file=Path(args.input_folder_slopes) / 'slopes_by_definition.tsv',
        sj_file=Path(args.input_folder_sj) / 'SJ.out.tab')
    selected_slopes_by_definition.to_csv(Path(args.output_folder) / 'selected_slopes_by_definition.tsv',
                                         sep='\t',
                                         index=False)

    selected_slopes_read_pairs = select_intron_slopes_by_sj_evidence(
        slopes_file=Path(args.input_folder_slopes) / 'slopes_read_pairs.tsv',
        sj_file=Path(args.input_folder_sj) / 'SJ.out.tab')
    selected_slopes_read_pairs.to_csv(Path(args.output_folder) / 'selected_slopes_read_pairs.tsv',
                                      sep='\t',
                                      index=False)

    selected_slopes_nascent_introns = select_intron_slopes_by_sj_evidence(
        slopes_file=Path(args.input_folder_slopes) / 'slopes_nascent_introns.tsv',
        sj_file=Path(args.input_folder_sj) / 'SJ.out.tab')
    selected_slopes_nascent_introns.to_csv(Path(args.output_folder) / 'selected_slopes_nascent_introns.tsv',
                                           sep='\t',
                                           index=False)

    selected_slopes_read_pairs_cummax = select_intron_slopes_by_sj_evidence(
        slopes_file=Path(args.input_folder_slopes) / 'slopes_read_pairs_cummax.tsv',
        sj_file=Path(args.input_folder_sj) / 'SJ.out.tab')
    selected_slopes_read_pairs_cummax.to_csv(Path(args.output_folder) / 'selected_slopes_read_pairs_cummax.tsv',
                                             sep='\t',
                                             index=False)
