import pandas as pd
from pathlib import Path
import argparse
import json


def aggregate_adapters(input_folder_path: Path, output_folder_path: Path) -> None:
    """
    Check adapters detected by Atria, and compare whether they are same across all samples.
    Saves the result to a JSON file.
    :param input_folder_path: Path to input folder; it should contain one subfolder per sample.
    :param output_folder_path: Folder to which the resulting JSON file will be written.
    :return: None.
    """
    adapters_1: set[str] = set()
    adapters_2: set[str] = set()
    for sub_folder in input_folder_path.iterdir():
        data = pd.read_csv(sub_folder / 'detected_adapters.tsv', delimiter='\t').squeeze()
        adapters_1.add(data['best_adapter1'])
        adapters_2.add(data['best_adapter2'])

    adapter_1 = adapters_1.pop() if len(adapters_1) == 1 else None
    adapter_2 = adapters_2.pop() if len(adapters_2) == 1 else None
    output_jon = {'adapter_1': adapter_1,
                  'adapter_2': adapter_2,
                  'adapters_successfully_detected': (adapter_1 is not None) and (adapter_2 is not None)}

    with open(output_folder_path / 'aggregated_adapters.json', 'w') as output_file:
        json.dump(output_jon, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        help='Folder containing output of Atria adapters detection.')
    parser.add_argument('--output_folder',
                        help='Folder to which the result will be saved.')
    args = parser.parse_args()
    aggregate_adapters(input_folder_path=Path(args.input_folder),
                       output_folder_path=Path(args.output_folder))
