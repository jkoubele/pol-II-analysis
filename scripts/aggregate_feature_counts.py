import argparse
from pathlib import Path
from typing import Optional
from tqdm import tqdm
import pandas as pd


def aggregate_feature_counts(input_folder_path: Path,
                             output_folder_path: Path,
                             count_column_name='/input_folder/Aligned.sortedByCoord.out.bam') -> None:
    """
    Aggregate individual feature count .tsv file into one.
    :param input_folder_path: Path to input folder; it should contain one subfolder per sample.
    :param output_folder_path: Folder to which the resulting .tsv file will be written.
    :return: None.
    """

    sample_to_counts: dict[str, pd.Series] = {}
    template_df: Optional[pd.DataFrame] = None
    for sample_folder in tqdm(input_folder_path.iterdir()):
        sample_df = pd.read_csv(sample_folder / 'feature_counts.tsv', delimiter='\t', header=1)
        sample_to_counts[sample_folder.name] = sample_df[count_column_name]
        if template_df is None:
            template_df = sample_df.drop(columns=count_column_name)

    for sample_name, counts in sample_to_counts.items():
        template_df[sample_name] = counts
    output_folder_path.mkdir(exist_ok=True, parents=True)
    template_df.to_csv(output_folder_path / 'aggregated_feature_counts.tsv', sep='\t', index=False)


if __name__ == "__main__":
    aggregate_feature_counts(input_folder_path=Path('/cellfile/datapublic/jkoubele/senescent_cells/feature_counts_exons'),
                             output_folder_path=Path(
                                 '/cellfile/datapublic/jkoubele/senescent_cells/aggregated_feature_counts_exons'))

    # parser = argparse.ArgumentParser()
    # parser.add_argument('--input_folder',
    #                     help='Folder containing count matrices produced by featureCounts .')
    # parser.add_argument('--output_folder',
    #                     help='Folder to which the result (single aggregated matrix) will be saved.')
    # args = parser.parse_args()
    # aggregate_feature_counts(input_folder_path=Path(args.input_folder),
    #                          output_folder_path=Path(args.output_folder))
    #
    #