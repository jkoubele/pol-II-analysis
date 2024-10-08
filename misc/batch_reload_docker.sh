#!/bin/bash

usage() {
    echo "Usage: $0 [-d <docker_image_path>] [-l <slurm_log_folder>]"
    exit 1
}

script_directory="$(cd "$(dirname "$0")" && pwd)"
repository_path="$(dirname "$script_directory")"

docker_image_path="$repository_path"/docker_images/bioinfo_tools.tar
slurm_log_folder="$repository_path"/slurm_logs

# Parse command line arguments
while getopts ":d:l:" opt; do
    case ${opt} in
        d )
            docker_image_path=$OPTARG
            ;;
        l )
            slurm_log_folder=$OPTARG
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

node_list=$(sinfo -N -h -o "%N" | sort | uniq)

for node in $node_list; do
    echo "Submitting job to reload docker on node: $node"
    sbatch --output="$slurm_log_folder"/%j_%x.log --error="$slurm_log_folder"/%j_%x.err \
    --nodelist="$node" "$repository_path"/misc/reload_docker_on_node.sh -d "$docker_image_path"
done
