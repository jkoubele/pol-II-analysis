usage() {
    echo "Usage: $0 -f <fai_file_name>"
    exit 1
}

# Variables to hold argument
fai_file_name=""

# Parse command line argument
while getopts ":f:" opt; do
    case ${opt} in
        f )
            fai_file_name=$OPTARG
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

# Check if mandatory argument is provided
if [ -z "$fai_file_name" ]; then
    echo "Error: Missing mandatory argument"
    usage
fi

for bed_file_name in forward_pairs reverse_pairs forward_nascent_introns reverse_nascent_introns;
do
  sort -k 1,1 -k 2,2n output_folder/$bed_file_name.bed > output_folder/${bed_file_name}_sorted.bed
  mv output_folder/${bed_file_name}_sorted.bed output_folder/$bed_file_name.bed
  bedtools genomecov -bga -split -i /output_folder/$bed_file_name.bed -g /genome_folder/"$fai_file_name" > \
  /output_folder/coverage_${bed_file_name}.bedGraph
  gzip output_folder/$bed_file_name.bed
done