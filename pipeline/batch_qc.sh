#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_folder> -o <output_folder> -d <docker_image_path> "
    echo "-l <log_folder> -e <error_folder>"
    exit 1
}

# Variables to hold arguments
input_folder=""
output_folder=""
docker_image_path=""

log_folder=""
error_folder=""

# Parse command line arguments
while getopts ":i:o:d:l:e:" opt; do
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
        l )
            log_folder=$OPTARG
            ;;
        e )
            error_folder=$OPTARG
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
if [ -z "$input_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] || \
[ -z "$log_folder" ] || [ -z "$error_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

for sub_folder in $input_folder; do
  file_name=$(basename "$sub_folder")
  dir_name=$(dirname "$sub_folder")
  echo "$file_name"
  echo "$dir_name"
done
