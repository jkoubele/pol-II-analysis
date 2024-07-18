#!/bin/bash

#SBATCH --job-name=feature_counts
#SBATCH --partition=all
#SBATCH --ntasks=12

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path>"
    echo " -g <genome_folder> -a <annotation_gtf_file_name> -s <strandedness>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""
genome_folder=""
annotation_gtf_file_name=""
strandedness=""

# Parse command line arguments
while getopts ":i:o:d:a:g:s:" opt; do
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
        g )
            genome_folder=$OPTARG
            ;;
        a )
            annotation_gtf_file_name=$OPTARG
            ;;
        s )
            strandedness=$OPTARG
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
|| [ -z "$genome_folder" ] || [ -z "$annotation_gtf_file_name" ] || [ -z "$strandedness" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run featureCounts
docker run --rm -v "$input_folder":/input_folder -v "$output_folder":/output_folder \
-v "$genome_folder":/genome_folder --security-opt seccomp=unconfined \
bioinfo_tools /bin/sh -c "featureCounts -p --countReadPairs -s $strandedness -T 12 \
-t gene -a /genome_folder/$annotation_gtf_file_name \
-o /output_folder/feature_counts.tsv \
/input_folder/Aligned.sortedByCoord.out.bam;
chmod 777 -R /output_folder"
