#!/bin/bash

#SBATCH --job-name=star_indexing
#SBATCH --ntasks=18

# Function to display usage information
usage() {
    echo "Usage: $0 -g <genome_folder> -a <annotation_gtf_file_name> -f <fasta_file_name>"
    echo "[-d <docker_image_path>]"
    exit 1
}

# Variables to hold arguments
genome_folder=""
annotation_gtf_file_name=""
fasta_file_name=""

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar

# Parse command line arguments
while getopts ":g:a:f:d:" opt; do
    case ${opt} in
        g )
            genome_folder=$OPTARG
            ;;
        a )
            annotation_gtf_file_name=$OPTARG
            ;;
        f )
            fasta_file_name=$OPTARG
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
if [ -z "$genome_folder" ] || [ -z "$annotation_gtf_file_name" ] || [ -z "$fasta_file_name" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Run STAR indexing
docker run --rm -v "$genome_folder":/genome_folder \
--security-opt seccomp=unconfined bioinfo_tools /bin/sh -c "STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir /genome_folder/STAR_index \
--genomeFastaFiles /genome_folder/$fasta_file_name
--sjdbGTFfile /genome_folder/$annotation_gtf_file_name
--outTmpDir /genome_folder/tmp \
--limitGenomeGenerateRAM 108000000000; \
chmod 777 -R /genome_folder"

