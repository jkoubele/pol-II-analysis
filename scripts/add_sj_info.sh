#!/bin/bash

#SBATCH --job-name=sj_info
#SBATCH --ntasks=15

# Function to display usage information
usage() {
    echo "Usage: $0 -i <intron_slopes_folder> -s <sj_folder> -o <output_folder>"
    echo " -d <docker_image_path> -c <script_folder> "
    exit 1
}

# Variables to hold arguments
intron_slopes_folder=""
sj_folder=""
output_folder=""
docker_image_path=""
script_folder=""



# Parse command line arguments
while getopts ":i:s:o:d:c:" opt; do
    case ${opt} in
        i )
            intron_slopes_folder=$OPTARG
            ;;
        o )
            output_folder=$OPTARG
            ;;
        d )
            docker_image_path=$OPTARG
            ;;
        s )
            sj_folder=$OPTARG
            ;;
        c )
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
if [ -z "$intron_slopes_folder" ] || [ -z "$output_folder" ] || [ -z "$docker_image_path" ] \
|| [ -z "$script_folder" ] || [ -z "$sj_folder" ]; then
    echo "Error: Missing mandatory arguments"
    usage
fi

# Check if the docker image is available, and load it from disk if it's not
if ! docker images --format "{{.Repository}}" | grep -q "^bioinfo_tools$"; then
    docker load -i "$docker_image_path"
fi

# Create output folder if it doesn't exist
mkdir "$output_folder" -p

# Run coverage computation
docker run --rm \
-v "$intron_slopes_folder":/intron_slopes_folder \
-v "$sj_folder":/sj_folder \
-v "$output_folder":/output_folder \
-v "$script_folder":/script_folder \
--security-opt seccomp=unconfined bioinfo_tools \
/bin/sh -c "python3 /script_folder/add_sj_info.py \
--input_folder_slopes /intron_slopes_folder \
--input_folder_sj /sj_folder \
--output_folder /output_folder; \
chmod 777 -R /output_folder"
