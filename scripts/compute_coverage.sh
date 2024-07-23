#!/bin/bash

#SBATCH --job-name=compute_coverage
#SBATCH --ntasks=12

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path>"
    echo " -s <strandedness> -c <script_folder> -g <genome_folder> -f <fai_file_name>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
strandedness=""
script_folder=""
genome_folder=""
fai_file_name=""


# Parse command line arguments
while getopts ":i:o:d:s:c:g:f:" opt; do
    case ${opt} in
        i )
            input_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        s )
            strandedness=$OPTARG
            ;;
        c )
            script_folder=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
            ;;
        f )
            fai_file_name=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] \
|| [ -z "$script_folder" ] || [ -z "$strandedness" ] || [ -z "$genome_folder" ] || [ -z "$fai_file_name" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Aggregate adapters
docker run --rm -v "$input_folder":/input_folder -v "$output_folder":/output_folder \
-v "$script_folder":/script_folder -v "$genome_folder":/genome_folder --security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "python3 /script_folder/extract_pairs_by_strand.py \
--input_folder /input_folder --output_folder /output_folder --strandendess_type $strandedness;  \
bedtools genomecov -bga -split -i /output_folder/forward.bed.gz -g /genome_folder/$fai_file_name > \
/output_folder/coverage_forward.bedGraph; \
bedtools genomecov -bga -split -i /output_folder/reverse.bed.gz -g /genome_folder/$fai_file_name > \
/output_folder/coverage_reverse.bedGraph; \
chmod 777 -R /output_folder"
