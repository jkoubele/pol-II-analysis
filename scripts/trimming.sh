#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --partition=all
#SBATCH --ntasks=15

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -a <aggregated_adapters_folder> -d <docker_image_path> "
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
aggregated_adapters_folder=""

# Parse command line arguments
while getopts ":i:o:d:a:" opt; do
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
        a )
            aggregated_adapters_folder=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] \
|| [ -z "$docker_image_path" ] || [ -z "$aggregated_adapters_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run Atria for trimming
docker run --rm -v "$input_folder":/input_folder -v "$output_folder":/output_folder \
-v "$aggregated_adapters_folder":/aggregated_adapters_folder --security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "atria --read1 /input_folder/R1.fastq.gz --read2 /input_folder/R2.fastq.gz \
--adapter1 \$(jq -r '.adapter_1' /aggregated_adapters_folder/aggregated_adapters.json) \
--adapter2 \$(jq -r '.adapter_2' /aggregated_adapters_folder/aggregated_adapters.json) \
--output-dir /output_folder -t 12; \
mv /output_folder/R1.atria.fastq.gz /output_folder/R1.fastq.gz; \
mv /output_folder/R2.atria.fastq.gz /output_folder/R2.fastq.gz; \
chmod 777 -R /output_folder"
