#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <intron_slopes_folder> -s <sj_folder> -o <output_folder>"
    echo "[-d <docker_image_path>] [-l <slurm_log_folder>] [-L]"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"

# Variables to hold arguments
intron_slopes_folder=""
sj_folder=""
output_folder=""

run_locally=false
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
slurm_log_folder="$repository_path"/slurm_logs


# Parse command line arguments
while getopts ":i:s:o:d:l:L" opt; do
    case ${opt} in
        i )
            intron_slopes_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        s )
            sj_folder=$OPTARG
            ;;
        l )
            slurm_log_folder=$OPTARG
            ;;
        L )
            run_locally=true
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            usage
            ;;
    esac
done

# Check if mandatory arguments are provided
if [ -z "$intron_slopes_folder" ] || [ -z "$output_folder" ] || [ -z "$sj_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

for sub_folder in "$intron_slopes_folder"/*; do
  sample_name=$(basename "$sub_folder")
  if [ "$run_locally" = true ]; then
    echo "Processing sample $sample_name"
    sh "$repository_path"/scripts/add_sj_info.sh -i "$sub_folder" -o "$output_folder"/"$sample_name" \
    -d "$docker_image_path" -s "$sj_folder"/"$sample_name" -c "$repository_path"/scripts
  else
    echo "Submitting sample $sample_name"
    sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
    "$repository_path"/scripts/add_sj_info.sh -i "$sub_folder" -o "$output_folder"/"$sample_name" \
    -d "$docker_image_path" -s "$sj_folder"/"$sample_name" -c "$repository_path"/scripts
  fi
done
