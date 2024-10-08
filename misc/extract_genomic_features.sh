#!/bin/bash

#SBATCH --job-name=extract_genomic_features
#SBATCH --ntasks=12

# Function to display usage information
usage() {
    echo "Usage: $0 -g <genome_folder> -a <annotation_gtf_file_name> -t <gtf_source_type> "
    echo "[-d <docker_image_path>] [-s <script_folder>]"
    exit 1
}

# Variables to hold arguments
genome_folder=""
annotation_gtf_file_name=""
gtf_source_type=""

script_folder="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_folder")"
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar

# Parse command line arguments
while getopts ":g:a:d:s:t:" opt; do
    case ${opt} in
        g )
            genome_folder=$OPTARG
            ;;
        a )
            annotation_gtf_file_name=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        s )
            script_folder=$OPTARG
            ;;
        t )
            gtf_source_type=$OPTARG
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
if [ -z "$genome_folder" ] || [ -z "$annotation_gtf_file_name" ] || [ -z "$gtf_source_type" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi


docker run --rm -v "$genome_folder":/genome_folder -v "$script_folder":/script_folder \
--security-opt seccomp=unconfined bioinfo_tools /bin/sh -c "python3 /script_folder/extract_genomic_features.py \
--genome_folder /genome_folder \
--gtf_file_name $annotation_gtf_file_name \
--gtf_source $gtf_source_type;  \
chmod 777 -R /genome_folder"
