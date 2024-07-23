#!/bin/bash

#SBATCH --job-name=extract_genomic_features
#SBATCH --ntasks=12

# Function to display usage information
usage() {
    echo "Usage: $0 -g <genome_folder> -a <annotation_gtf_file_name> -d <docker_image_path> -s <script_folder>"
    exit 1
}

# Variables to hold arguments
genome_folder=""
annotation_gtf_file_name=""
docker_image_path=""
script_folder=""

# Parse command line arguments
while getopts ":g:a:d:s:" opt; do
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
if [ -z "$genome_folder" ] || [ -z "$annotation_gtf_file_name" ] || [ -z "$docker_image_path" ] \
|| [ -z "$script_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi


docker run --rm -v "$genome_folder":/genome_folder -v "$script_folder":/script_folder \
--security-opt seccomp=unconfined bioinfo_tools /bin/sh -c "python3 /script_folder/extract_genomic_features.py \
--genome_folder /genome_folder --gtf_file_name $annotation_gtf_file_name;  \
chmod 777 -R /genome_folder"
