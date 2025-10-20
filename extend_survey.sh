#!/bin/bash

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to display help message
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Description:
    Process duplicated gene pairs by analyzing flanking sequences using promer alignment.
    This script extracts upstream/downstream sequences and performs comparative analysis.

Required Options:
    -fd, --first-dir <dir>          First directory name (e.g., 5k)
    -sd, --second-dir <dir>         Second directory name (e.g., 10k)
    -f,  --flank-size <size>        Flanking sequence size in bp (e.g., 10000)
    -pt, --promer-threads <num>     Number of threads for promer (e.g., 5)
    -gt, --general-threads <num>    Number of threads for general processing (e.g., 20)
    -it, --identity <threshold>     Identity threshold for delta-filter (e.g., 90)
    -d,  --distance <dist>          Distance parameter for up/down processing (e.g., 5000)

Optional:
    -h,  --help                     Display this help message and exit

Environment Variables Required:
    \$workspace    Path to workspace directory
    \$genome       Path to genome FASTA file
    \$script       Path to script directory (containing up.sh and down.sh)

Example:
    export workspace=/path/to/workspace
    export genome=/path/to/genome.fa
    export script=/path/to/scripts
    
    $(basename "$0") -fd 5k -sd 10k -f 10000 -pt 5 -gt 20 -it 90 -d 5000

EOF
    exit 0
}

# Function to print error messages
print_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Function to print info messages
print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

# Function to print warning messages
print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Initialize variables
first_dir=""
second_dir=""
flank_size=""
cpu1=""
cpu2=""
identity_threshold=""
distance=""

# Parse command line arguments
if [ $# -eq 0 ]; then
    show_help
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            ;;
        -fd|--first-dir)
            first_dir="$2"
            shift 2
            ;;
        -sd|--second-dir)
            second_dir="$2"
            shift 2
            ;;
        -f|--flank-size)
            flank_size="$2"
            shift 2
            ;;
        -pt|--promer-threads)
            cpu1="$2"
            shift 2
            ;;
        -gt|--general-threads)
            cpu2="$2"
            shift 2
            ;;
        -it|--identity)
            identity_threshold="$2"
            shift 2
            ;;
        -d|--distance)
            distance="$2"
            shift 2
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$first_dir" ] || [ -z "$second_dir" ] || [ -z "$flank_size" ] || \
   [ -z "$cpu1" ] || [ -z "$cpu2" ] || [ -z "$identity_threshold" ] || [ -z "$distance" ]; then
    print_error "Missing required parameters"
    echo "Use -h or --help for usage information"
    exit 1
fi

# Validate environment variables
if [ -z "$workspace" ]; then
    print_error "Environment variable \$workspace is not set"
    exit 1
fi

if [ -z "$genome" ]; then
    print_error "Environment variable \$genome is not set"
    exit 1
fi

if [ -z "$script" ]; then
    print_error "Environment variable \$script is not set"
    exit 1
fi

# Check if required files exist
if [ ! -f "$genome" ]; then
    print_error "Genome file not found: $genome"
    exit 1
fi

if [ ! -d "$workspace/Dupgene/process/${first_dir}" ]; then
    print_error "First directory not found: $workspace/Dupgene/process/${first_dir}"
    exit 1
fi

# Print configuration
print_info "Starting duplicated gene processing pipeline"
print_info "Configuration:"
echo "  First directory:      $first_dir"
echo "  Second directory:     $second_dir"
echo "  Flank size:           $flank_size bp"
echo "  Promer threads:       $cpu1"
echo "  General threads:      $cpu2"
echo "  Identity threshold:   $identity_threshold"
echo "  Distance:             $distance"
echo "  Workspace:            $workspace"
echo ""

# Main processing pipeline
print_info "Step 1: Preparing directories"
cd "$workspace/Dupgene/process/${first_dir}" || exit 1
rm -rf "$workspace/Dupgene/process/${second_dir}"
mkdir -p "$workspace/Dupgene/process/${second_dir}"

print_info "Step 2: Extracting gene IDs for upstream and downstream analysis"
awk '$3=="up"{print $4}' "${first_dir}_up_down_check" | sort -u > "$workspace/Dupgene/process/${second_dir}/${second_dir}_up_list"
awk '$3=="down"{print $4}' "${first_dir}_up_down_check" | sort -u > "$workspace/Dupgene/process/${second_dir}/${second_dir}_down_list"

cd "$workspace/Dupgene/process/${second_dir}" || exit 1
cat "${second_dir}_up_list" "${second_dir}_down_list" | sort -u > "${second_dir}_list"

print_info "Step 3: Extracting flanking sequences"
# Extract upstream sequences
grep -wf "${second_dir}_list" "$workspace/DupGen_finder/Dupgene.pairs" | awk '{print $2"\n"$3}' | sort -u | \
grep -wf - "$workspace/DupGen_finder/Dupgene.bed" | \
bedtools flank -i - -g "$workspace/Dupgene/genome.fa.len" -l "$flank_size" -r 0 -s > up.Dupgene.bed

# Extract downstream sequences
grep -wf "${second_dir}_list" "$workspace/DupGen_finder/Dupgene.pairs" | awk '{print $2"\n"$3}' | sort -u | \
grep -wf - "$workspace/DupGen_finder/Dupgene.bed" | \
bedtools flank -i - -g "$workspace/Dupgene/genome.fa.len" -l 0 -r "$flank_size" -s > down.Dupgene.bed

print_info "Step 4: Extracting FASTA sequences from genome"
bedtools getfasta -s -fi "$genome" -bed up.Dupgene.bed -fo up.Dupgene.temp.fa -name
bedtools getfasta -s -fi "$genome" -bed down.Dupgene.bed -fo down.Dupgene.temp.fa -name

# Clean FASTA headers
sed 's/::/\t/' up.Dupgene.temp.fa | cut -f 1 > up.Dupgene.fa
sed 's/::/\t/' down.Dupgene.temp.fa | cut -f 1 > down.Dupgene.fa
rm up.Dupgene.temp.fa down.Dupgene.temp.fa

print_info "Step 5: Preparing temporary directory and generating gene ID files"
rm -rf "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP"
mkdir -p "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP"
cd "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP" || exit 1

rm -f get_id.sh* 2>/dev/null
while read -r i; do
    echo "mkdir -p $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
cd $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
awk -va=\"$i\" '\$1==a {print \$2\"\\n\"\$3}' $workspace/DupGen_finder/Dupgene.pairs > ${i}.id"
done < "../${second_dir}_list" > get_id.sh
ParaFly -c get_id.sh -CPU "${cpu2}"

print_info "Step 6: Running promer analysis for upstream sequences"
rm -f promer.up.sh* 2>/dev/null
while read -r i; do
    echo "cd $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
sed -n '1p' ${i}.id > a.id; \
sed -n '2p' ${i}.id > b.id; \
seqkit grep -f a.id ../../up.Dupgene.fa > a_up_${second_dir}.fa; \
seqkit grep -f b.id ../../up.Dupgene.fa > b_up_${second_dir}.fa; \
\$MUMmer/promer --maxmatch a_up_${second_dir}.fa b_up_${second_dir}.fa -p up; \
\$MUMmer/delta-filter -i ${identity_threshold} up.delta > up_filtered.delta; \
\$MUMmer/show-coords -Trcdl up_filtered.delta > up_filtered.coords; \
\$MUMmer/mummerplot --png up_filtered.delta -p up"
done < "../${second_dir}_up_list" > promer.up.sh
ParaFly -c promer.up.sh -CPU "${cpu1}" > up_script.log 2>&1

print_info "Step 7: Running promer analysis for downstream sequences"
rm -f promer.down.sh* 2>/dev/null
while read -r i; do
    echo "cd $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
sed -n '1p' ${i}.id > a.id; \
sed -n '2p' ${i}.id > b.id; \
seqkit grep -f a.id ../../down.Dupgene.fa > a_down_${second_dir}.fa; \
seqkit grep -f b.id ../../down.Dupgene.fa > b_down_${second_dir}.fa; \
\$MUMmer/promer --maxmatch a_down_${second_dir}.fa b_down_${second_dir}.fa -p down; \
\$MUMmer/delta-filter -i ${identity_threshold} down.delta > down_filtered.delta; \
\$MUMmer/show-coords -Trcdl down_filtered.delta > down_filtered.coords; \
\$MUMmer/mummerplot --png down_filtered.delta -p down"
done < "../${second_dir}_down_list" > promer.down.sh
ParaFly -c promer.down.sh -CPU "${cpu1}" > down_script.log 2>&1

print_info "Step 8: Processing upstream results"
cd "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP" || exit 1
rm -f run_up.sh* 2>/dev/null
while read -r i; do
    echo "cd $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
export workspace=$workspace; \
export script=$script; \
\$script/up.sh -p ${i} -d ${distance} -l ${flank_size} 2>/dev/null"
done < <(cut -f 1 "../${second_dir}_up_list") > run_up.sh
ParaFly -c run_up.sh -CPU "${cpu2}"

print_info "Step 9: Processing downstream results"
cd "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP" || exit 1
rm -f run_down.sh* 2>/dev/null
while read -r i; do
    echo "cd $workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP/${i}; \
export workspace=$workspace; \
export script=$script; \
\$script/down.sh -p ${i} -d ${distance} -l ${flank_size} 2>/dev/null"
done < <(cut -f 1 "../${second_dir}_down_list") > run_down.sh
ParaFly -c run_down.sh -CPU "${cpu2}"

print_info "Step 10: Collecting results"
cd "$workspace/Dupgene/process/${second_dir}/${second_dir}_TEMP" || exit 1
find . -name "check_up" -exec cat {} + | awk 'NF==4' > "../${second_dir}_up_check"
find . -name "check_down" -exec cat {} + | awk 'NF==4' > "../${second_dir}_down_check"
cat "../${second_dir}_up_check" "../${second_dir}_down_check" > "../${second_dir}_up_down_check"

print_info "Pipeline completed successfully!"
print_info "Results saved to: $workspace/Dupgene/process/${second_dir}/${second_dir}_up_down_check"