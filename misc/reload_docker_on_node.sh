#!/bin/bash

#SBATCH --job-name=reload_docker

usage() {
    echo "Usage: $0 [-d <docker_image_path>]"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"
docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar

# Parse command line arguments
while getopts ":d:" opt; do
    case ${opt} in
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

docker image prune -a -f
docker load -i "$docker_image_path"
