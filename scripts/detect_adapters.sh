#!/bin/bash

#SBATCH --job-name=detect_adapters
#SBATCH --partition=all
#SBATCH --ntasks=3

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""

# Parse command line arguments
while getopts ":i:o:d:" opt; do
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run Atria to detect adapters
docker run --rm -v "$input_folder":/input_folder -v "$output_folder":/output_folder \
--security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "atria --read1 /input_folder/R1.fastq.gz --read2 /input_folder/R2.fastq.gz --detect-adapter; \
mv atria_adapter_detect_summary* /output_folder/detected_adapters.tsv; \
chmod 777 -R /output_folder"
