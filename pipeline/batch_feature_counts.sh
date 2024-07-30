#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -g <genome_folder>"
    echo " -a <annotation_gtf_file_name> -s <strandedness> [-f <feature_type>]"
    echo "[-d <docker_image_path>] [-l <slurm_log_folder>] [-L]"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"

# Variables to hold arguments
input_folder=""
output_folder=""
genome_folder=""
annotation_gtf_file_name=""
strandedness=""
feature_type="gene"
run_locally=false
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
slurm_log_folder="$repository_path"/slurm_logs

# Parse command line arguments
while getopts ":i:o:g:d:a:s:f:l:L" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        a )
            annotation_gtf_file_name=$OPTARG
            ;;
        s )
            strandedness=$OPTARG
            ;;
        f )
            feature_type=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$genome_folder" ] \
|| [ -z "$annotation_gtf_file_name" ] || [ -z "$strandedness" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

for sub_folder in "$input_folder"/*; do
  sample_name=$(basename "$sub_folder")
  if [ "$run_locally" = true ]; then
    echo "Processing sample $sample_name"
    sh "$repository_path"/scripts/feature_counts.sh -i "$sub_folder" -o "$output_folder"/"$sample_name" \
    -d "$docker_image_path" -g "$genome_folder" -a "$annotation_gtf_file_name" -s "$strandedness" -f "$feature_type"
  else
    echo "Submitting sample $sample_name"
    sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
    "$repository_path"/scripts/feature_counts.sh -i "$sub_folder" -o "$output_folder"/"$sample_name" \
    -d "$docker_image_path" -g "$genome_folder" -a "$annotation_gtf_file_name" -s "$strandedness" -f "$feature_type"
  fi
done
