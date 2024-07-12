# Function to display help
show_help() {
  echo "Usage: $0 [-o output_folder] [--help]"
  echo ""
  echo "Options:"
  echo "  -o      Specify the output folder (default is '../docker_images')."
  echo "  --help  Display this help message."
}

# Check for --help before parsing other options
for arg in "$@"; do
  case $arg in
    --help)
      show_help
      exit 0
      ;;
  esac
done

# default value of output_folder
output_folder="$(dirname "$0")/docker_images"

# Parse command line options
while getopts ":o:" opt; do
  case ${opt} in
    o )
      output_folder=$OPTARG
      ;;
    \? )
      echo "Usage: cmd [-o output_folder]"
      exit 1
      ;;
  esac
done

# Script logic below
mkdir "$output_folder" -p

docker build -t bioinfo_tools "$(dirname "$0")/dockerfiles/bioinfo_tools"
docker save -o "$output_folder"/bioinfo_tools.tar bioinfo_tools
