#!/bin/bash

#SBATCH --job-name=estimate_intron_slopes
#SBATCH --ntasks=15

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path>"
    echo " -s <script_folder> -g <genome_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
script_folder=""
genome_folder=""

# Parse command line arguments
while getopts ":i:o:d:s:g:" opt; do
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
            script_folder=$OPTARG
            ;;
        g )
            genome_folder=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ]||
[ -z "$script_folder" ] || [ -z "$genome_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_r$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run estimation of slopes in docker
docker run --rm -v "$input_folder":/input_folder \
-v "$output_folder":/output_folder \
-v "$script_folder":/script_folder \
-v "$genome_folder":/genome_folder \
--security-opt seccomp=unconfined \
bioinfo_r /bin/sh -c "Rscript /script_folder/estimate_intron_slopes.R \
--input_folder /input_folder \
--output_folder /output_folder \
--introns_file /genome_folder/introns.bed; \
chmod 777 -R /output_folder"
