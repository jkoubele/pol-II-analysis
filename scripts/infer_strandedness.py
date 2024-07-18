import argparse
import json
from pathlib import Path

import pandas as pd


def infer_strandedness(input_folder: Path, output_folder: Path) -> None:
    """
    Infer strandedness of samples from feature counts files produced by STAR and
    writes output to JSON.
    :param input_folder: Input folder, containing subfolder for each sample. Each subfolder
    is expected to be the result of STAR alignment and contains the file 'ReadsPerGene.out.tab'.
    :param output_folder: Output folder to which output JSON will be written
    :return: None.
    """
    num_stranded_1 = 0
    num_stranded_2 = 0
    for subfolder in input_folder.iterdir():
        gene_counts = pd.read_csv(subfolder / 'ReadsPerGene.out.tab', sep='\t',
                                  names=['feature', 'unstranded', 'stranded', 'reversed'])
        gene_counts = gene_counts.set_index('feature')
        if gene_counts.loc['N_unmapped']['stranded'] < gene_counts.loc['N_unmapped']['reversed']:
            num_stranded_1 += 1
        else:
            num_stranded_2 += 1

    strandendess_type = '0'
    if num_stranded_1 == 0:
        strandendess_type = '2'
    elif num_stranded_2 == 0:
        strandendess_type = '1'
    output_json = {'strandendess_type': strandendess_type,
                   'num_stranded_1': num_stranded_1,
                   'num_stranded_2': num_stranded_2}

    with open(output_folder / 'strandedness_info.json', 'w') as output_file:
        json.dump(output_json, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        help='Folder containing subfolders with output of STAR alignment.')
    parser.add_argument('--output_folder',
                        help='Folder to which the result will be saved.')
    args = parser.parse_args()
    infer_strandedness(input_folder=Path(args.input_folder),
                       output_folder=Path(args.output_folder))
