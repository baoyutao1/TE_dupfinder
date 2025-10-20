#!/bin/bash

###################################################################################
# Script: calculate_kaks.sh
# Description: Calculate Ka/Ks ratios for duplicate gene pairs using ParaAT and 
#              KaKs_Calculator
###################################################################################

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to display usage information
usage() {
    cat << EOF
Usage: $(basename "$0") -w WORKSPACE -c CDS_FILE -p PEP_FILE [OPTIONS]

Required Arguments:
  -w, --workspace PATH    Workspace directory (must contain DupGen_finder/Dupgene.pairs)
  -c, --cds PATH          CDS sequence file in FASTA format
  -p, --pep PATH          Protein sequence file in FASTA format

Optional Arguments:
  -t, --threads NUM       Number of threads for parallel processing (default: 30)
  -h, --help              Display this help message and exit

Description:
  This script calculates Ka/Ks ratios for duplicate gene pairs identified by DupGen_finder.
  It performs the following steps:
    1. Prepares input files by cleaning gene IDs
    2. Aligns sequences using MAFFT via ParaAT
    3. Calculates Ka/Ks values using KaKs_Calculator (YN model)
    4. Filters and formats results (Ka<3, Ks<10)

Requirements:
  - ParaAT.pl (for sequence alignment)
  - MAFFT (multiple sequence alignment)
  - KaKs_Calculator (for Ka/Ks calculation)
  - DupGen_finder output (Dupgene.pairs file)

Output:
  - sub_dupgene.kaks.txt: Final Ka/Ks results with duplicate type classification

Example:
  $(basename "$0") -w /path/to/workspace -c genome.cds.fa -p genome.pep.faa -t 40

EOF
    exit 0
}

# Function to print error messages
error_exit() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

# Function to print info messages
info_msg() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

# Function to print warning messages
warn_msg() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

# Function to check if required commands exist
check_dependencies() {
    local deps=("ParaAT.pl" "mafft" "KaKs_Calculator")
    local missing=()
    
    for cmd in "${deps[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [ ${#missing[@]} -ne 0 ]; then
        error_exit "Missing required dependencies: ${missing[*]}"
    fi
}

# Function to validate input files
validate_inputs() {
    # Check workspace directory
    if [ ! -d "$WORKSPACE" ]; then
        error_exit "Workspace directory does not exist: $WORKSPACE"
    fi
    
    # Check for DupGen_finder output
    if [ ! -f "$WORKSPACE/DupGen_finder/Dupgene.pairs" ]; then
        error_exit "DupGen_finder output not found: $WORKSPACE/DupGen_finder/Dupgene.pairs"
    fi
    
    # Check CDS file
    if [ ! -f "$CDS_FILE" ]; then
        error_exit "CDS file not found: $CDS_FILE"
    fi
    
    # Check PEP file
    if [ ! -f "$PEP_FILE" ]; then
        error_exit "Protein file not found: $PEP_FILE"
    fi
}

# Initialize variables
WORKSPACE=""
CDS_FILE=""
PEP_FILE=""
THREADS=30

# Parse command line arguments
if [ $# -eq 0 ]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--workspace)
            WORKSPACE="$2"
            shift 2
            ;;
        -c|--cds)
            CDS_FILE="$2"
            shift 2
            ;;
        -p|--pep)
            PEP_FILE="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            error_exit "Unknown option: $1. Use -h or --help for usage information."
            ;;
    esac
done

# Validate required arguments
if [ -z "$WORKSPACE" ] || [ -z "$CDS_FILE" ] || [ -z "$PEP_FILE" ]; then
    error_exit "Missing required arguments. Use -h or --help for usage information."
fi

# Main execution starts here
info_msg "Starting Ka/Ks calculation pipeline..."
info_msg "Workspace: $WORKSPACE"
info_msg "CDS file: $CDS_FILE"
info_msg "Protein file: $PEP_FILE"
info_msg "Threads: $THREADS"

# Check dependencies
info_msg "Checking dependencies..."
check_dependencies

# Validate input files
info_msg "Validating input files..."
validate_inputs

# Create and navigate to kaks directory
KAKS_DIR="$WORKSPACE/kaks"
mkdir -p "$KAKS_DIR"
cd "$KAKS_DIR" || error_exit "Failed to change directory to $KAKS_DIR"

info_msg "Working directory: $KAKS_DIR"

# Step 1: Prepare input files
info_msg "Step 1/5: Preparing input files..."

# Clean gene IDs in duplicate pairs file (remove special characters)
cut -f 2,3 "$WORKSPACE/DupGen_finder/Dupgene.pairs" | \
    sed 's/\./_/g' | \
    sed -e 's/-/_/g' -e 's/|/_/g' -e 's/:/_/g' > input_dupgene.pairs || \
    error_exit "Failed to create input_dupgene.pairs"

# Create classification file with gene pairs
awk '{print $4"\t"$2"\t"$3}' "$WORKSPACE/DupGen_finder/Dupgene.pairs" | \
    sed 's/\./_/g' | \
    sed -e 's/-/_/g' -e 's/|/_/g' -e 's/:/_/g' | \
    awk '{print $1"\t"$2"-"$3}' > sub_dupgene.pairs || \
    error_exit "Failed to create sub_dupgene.pairs"

# Clean gene IDs in CDS file
sed '/>/ s/\./_/g' "$CDS_FILE" | \
    sed -e '/>/ s/-/_/g' -e '/>/ s/|/_/g' -e '/>/ s/:/_/g' > cds.fa || \
    error_exit "Failed to create cleaned CDS file"

# Clean gene IDs in protein file
sed '/>/ s/\./_/g' "$PEP_FILE" | \
    sed -e '/>/ s/-/_/g' -e '/>/ s/|/_/g' -e '/>/ s/:/_/g' > pep.faa || \
    error_exit "Failed to create cleaned protein file"

info_msg "Input files prepared successfully"

# Step 2: Run ParaAT for sequence alignment
info_msg "Step 2/5: Running ParaAT for sequence alignment..."

# Create processor configuration file
echo "$THREADS" > proc

# Remove old temporary directory if exists
rm -rf ./temp

# Run ParaAT: align homologous sequences and back-translate to CDS
$ParaAT/ParaAT.pl -h input_dupgene.pairs -n cds.fa -a pep.faa -p proc \
    -m mafft -f axt -g -k -o ./temp || \
    error_exit "ParaAT.pl failed"

info_msg "Sequence alignment completed"

# Step 3: Merge alignment results
info_msg "Step 3/5: Merging alignment results..."

find ./temp -name "*.axt" -exec cat {} \; > merge_align.axt || \
    error_exit "Failed to merge alignment files"

# Check if merged file is not empty
if [ ! -s merge_align.axt ]; then
    error_exit "No alignment results found"
fi

info_msg "Alignment files merged successfully"

# Step 4: Calculate Ka/Ks values
info_msg "Step 4/5: Calculating Ka/Ks values using YN model..."

$KaKs_Calculator/KaKs_Calculator -m YN -i merge_align.axt -o Calculator_result.txt || \
    error_exit "KaKs_Calculator failed"

info_msg "Ka/Ks calculation completed"

# Step 5: Filter and format results
info_msg "Step 5/5: Filtering and formatting results..."

# Extract columns, filter by Ka<3 and Ks<10, join with duplicate type, add header
cut -f 1,3,4,5 Calculator_result.txt | \
    awk '$3<3&&$4<10' | \
    awk 'NR==FNR{a[$1]=$0;next}{print $1"\t"a[$2]}' - sub_dupgene.pairs | \
    awk 'NF==5{print $0}' | \
    sed '1i Class\tSequence\tKa\tKs\tKa/Ks' > sub_dupgene.kaks.txt || \
    error_exit "Failed to create final output file"

# Count results
RESULT_COUNT=$(tail -n +2 sub_dupgene.kaks.txt | wc -l)

info_msg "Results filtered and formatted successfully"
info_msg "Total gene pairs with Ka/Ks values: $RESULT_COUNT"

# Summary
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Pipeline completed successfully!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Output files:"
echo "  - Final results: $KAKS_DIR/sub_dupgene.kaks.txt"
echo "  - Raw Ka/Ks: $KAKS_DIR/Calculator_result.txt"
echo "  - Merged alignments: $KAKS_DIR/merge_align.axt"
echo ""
echo "Filtering criteria applied:"
echo "  - Ka < 3"
echo "  - Ks < 10"
echo ""

exit 0