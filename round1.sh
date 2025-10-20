#!/bin/bash

################################################################################

# Display help information
show_help() {
    cat << EOF
Usage: $(basename "$0") -p <pair> [-d <distance>] [-l <length_threshold>] [-s <script_dir>] [-h]

Description:
    Process upstream and downstream genomic regions for gene pairs.
    Analyzes syntenic regions from coordinate files and filters results
    based on distance and length thresholds.

Required Options:
    -p <pair>              Gene pair identifier (required)

Optional Parameters:
    -d <distance>          Distance threshold for region filtering (default: 500)
    -l <length_threshold>  Length threshold for coverage checking (default: none)
    -s <script_dir>        Directory containing required Python scripts
                          (default: uses \$script environment variable)
    -h, --help            Display this help message and exit

Input Files (must exist in current directory):
    - up_filtered.coords     Upstream coordinate file
    - down_filtered.coords   Downstream coordinate file
    - ../../up.Dupgene.bed   Upstream gene annotation
    - ../../down.Dupgene.bed Downstream gene annotation

Output Files:
    - up.txt                 Merged upstream regions
    - down.txt               Merged downstream regions
    - check_up               Upstream coverage validation
    - check_down             Downstream coverage validation
    - Various intermediate .bed and .temp files

Required Python Scripts (in script directory):
    - filter_raw_data.py
    - filter_bed.py
    - up.py
    - down.py
    - merge_bed_regions.py

Examples:
    # Basic usage with required pair parameter
    $(basename "$0") -p gene_pair_01

    # Specify distance threshold
    $(basename "$0") -p gene_pair_01 -d 1000

    # Full parameter specification
    $(basename "$0") -p gene_pair_01 -d 1000 -l 5000 -s /path/to/scripts

Notes:
    - The script must be run from a directory containing the input coordinate files
    - Ensure all Python scripts are available in the specified script directory
    - Python 2 is required for the helper scripts

EOF
}

# Initialize variables with default values
distance=500
pair=""
length_threshold=""
script_dir="${script:-}"  # Use $script environment variable if set

# Parse command-line options
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_help
            exit 0
            ;;
        -d)
            distance="$2"
            shift 2
            ;;
        -p)
            pair="$2"
            shift 2
            ;;
        -l)
            length_threshold="$2"
            shift 2
            ;;
        -s)
            script_dir="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            echo "Use -h or --help for usage information" >&2
            exit 1
            ;;
    esac
done

# Validate required parameters
if [ -z "$pair" ]; then
    echo "Error: Missing required parameter -p <pair>" >&2
    echo "Use -h or --help for usage information" >&2
    exit 1
fi

# Validate script directory
if [ -z "$script_dir" ]; then
    echo "Error: Script directory not specified" >&2
    echo "Please set \$script environment variable or use -s option" >&2
    exit 1
fi

if [ ! -d "$script_dir" ]; then
    echo "Error: Script directory does not exist: $script_dir" >&2
    exit 1
fi

# Check for required input files
required_files=("up_filtered.coords" "down_filtered.coords" "../../up.Dupgene.bed" "../../down.Dupgene.bed")
for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "Error: Required input file not found: $file" >&2
        exit 1
    fi
done

# Set base directory
base_dir=$(pwd)
cd "$base_dir" || exit 1

echo "=========================================="
echo "Processing Gene Pair: $pair"
echo "Distance Threshold: $distance"
echo "Length Threshold: ${length_threshold:-not specified}"
echo "Script Directory: $script_dir"
echo "=========================================="

################################################################################
# UPSTREAM REGION PROCESSING
################################################################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting upstream processing..."

# Extract and format upstream coordinates
# Remove header lines (1-4) and format as BED: chr, start, end, gene_id
sed '1,4d' up_filtered.coords | awk '{print "chr\t"$1"\t"$2"\t"$16}' | \
    awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 1

sed '1,4d' up_filtered.coords | awk '{print "chr\t"$3"\t"$4"\t"$17}' | \
    awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 2

# Extract unique gene identifiers
cut -f 4 1 | sort -u > sp1
cut -f 4 2 | sort -u > sp2

# Calculate gene lengths from annotation files
length1=$(grep -wf sp1 ../../up.Dupgene.bed | awk '{print $3-$2}')
length2=$(grep -wf sp2 ../../up.Dupgene.bed | awk '{print $3-$2}')

echo "  Reference gene length: $length1"
echo "  Query gene length: $length2"

# Combine reference and query coordinates, sort by different criteria
paste 1 2 | sort -k 7nr -k 3nr | \
    awk -va=$length1 -vb=$length2 -vc=$distance \
    'NR==1 && (a-$3)<=c && (b-$7)<=c {flag=1} flag' > reference_query_up_raw.bed

paste 1 2 | sort -k 3nr -k 7nr | \
    awk -va=$length1 -vb=$length2 -vc=$distance \
    'NR==1 && (a-$3)<=c && (b-$7)<=c {flag=1} flag' > reference_query_up_raw_sort.bed

# Filter raw data to remove outliers (keep values within -500 to +500 relative range)
python2 "$script_dir/filter_raw_data.py" reference_query_up_raw.bed -d -500 -u 500 > temp1
python2 "$script_dir/filter_raw_data.py" reference_query_up_raw_sort.bed -d -500 -u 500 > temp2

# Merge filtered results and remove duplicates
cat temp1 temp2 | sort -u | sort -k 7nr -k 3nr > reference_query_up_clean.bed

# Split into separate reference and query files
cut -f 1-4 reference_query_up_clean.bed | sort -k 3nr -k 2nr > reference_up_raw.bed
cut -f 5-8 reference_query_up_clean.bed | sort -k 3nr -k 2nr > query_up_raw.bed

# Clean up temporary files
rm 1 2 temp1 temp2 reference_query_up_raw_sort.bed

# Remove completely contained regions (nested intervals)
python2 "$script_dir/filter_bed.py" reference_up_raw.bed > reference_up_clean.bed
python2 "$script_dir/filter_bed.py" query_up_raw.bed > query_up_clean.bed

# Process upstream regions with distance constraints
python2 "$script_dir/up.py" -c "$length1" -d "$distance" reference_up_clean.bed > reference_up.bed.temp
python2 "$script_dir/up.py" -c "$length2" -d "$distance" query_up_clean.bed > query_up.bed.temp

# Merge reference and query results, keeping only complete pairs
awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_up.bed.temp reference_query_up_raw.bed | \
    awk 'NF==9' | \
    awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_up.bed.temp - | \
    awk 'NF==10' | sort -k 3nr | cut -f1-4 > reference_up.temp

awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_up.bed.temp reference_query_up_raw.bed | \
    awk 'NF==9' | \
    awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_up.bed.temp - | \
    awk 'NF==10' | sort -k 7nr | cut -f5-8 > query_up.temp

# Merge overlapping regions and create final upstream output
python2 "$script_dir/merge_bed_regions.py" reference_up.temp | paste sp1 - | cut -f1,3-4 > sp1.up.temp
python2 "$script_dir/merge_bed_regions.py" query_up.temp | paste sp2 - | cut -f1,3-4 | cat sp1.up.temp - > up.txt

rm sp1.up.temp

# Validate upstream coverage against length threshold (90% coverage required)
if [ -n "$length_threshold" ]; then
    awk -v a="$length1" -v b="$length2" -v len_threshold="$length_threshold" \
        'a == len_threshold && b == len_threshold' up.txt | \
        awk '{print $1"\t"$3-$2+1"\tup"}' | \
        awk -va="$pair" -v len_threshold="$length_threshold" \
        '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_up
else
    # If no length threshold specified, create empty check file
    touch check_up
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Upstream processing completed"

################################################################################
# DOWNSTREAM REGION PROCESSING
################################################################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting downstream processing..."

# Extract and format downstream coordinates
sed '1,4d' down_filtered.coords | awk '{print "chr\t"$1"\t"$2"\t"$16}' | \
    awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 1

sed '1,4d' down_filtered.coords | awk '{print "chr\t"$3"\t"$4"\t"$17}' | \
    awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 2

# Combine and sort by start positions (ascending)
paste 1 2 | sort -k 6n -k 2n | \
    awk -va=$distance 'NR==1 && ($2<=a && $6<=a) {flag=1} flag' > reference_query_down_raw.bed

paste 1 2 | sort -k 2n -k 6n | \
    awk -va=$distance 'NR==1 && ($2<=a && $6<=a) {flag=1} flag' > reference_query_down_raw_sort.bed

# Extract unique gene identifiers for downstream
cut -f 4 1 | sort -u > sp1
cut -f 4 2 | sort -u > sp2

# Calculate downstream gene lengths
length3=$(grep -wf sp1 ../../down.Dupgene.bed | awk '{print $3-$2}')
length4=$(grep -wf sp2 ../../down.Dupgene.bed | awk '{print $3-$2}')

echo "  Reference gene length: $length3"
echo "  Query gene length: $length4"

# Filter downstream raw data
python2 "$script_dir/filter_raw_data.py" reference_query_down_raw.bed -d -500 -u 500 > temp1
python2 "$script_dir/filter_raw_data.py" reference_query_down_raw_sort.bed -d -500 -u 500 > temp2

# Merge and clean downstream data
cat temp1 temp2 | sort -u | sort -k 6n -k 2n > reference_query_down_clean.bed

# Split into reference and query files
cut -f 1-4 reference_query_down_clean.bed | sort -k 2n -k 3n > reference_down_raw.bed
cut -f 5-8 reference_query_down_clean.bed | sort -k 2n -k 3n > query_down_raw.bed

rm 1 2 temp1 temp2 reference_query_down_raw_sort.bed

# Filter nested regions for downstream
python2 "$script_dir/filter_bed.py" reference_down_raw.bed > reference_down_clean.bed
python2 "$script_dir/filter_bed.py" query_down_raw.bed > query_down_clean.bed

# Process downstream regions with distance threshold
python2 "$script_dir/down.py" reference_down_clean.bed -d "$distance" > reference_down.bed.temp
python2 "$script_dir/down.py" query_down_clean.bed -d "$distance" > query_down.bed.temp

# Merge downstream reference and query results
awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_down.bed.temp reference_query_down_raw.bed | \
    awk 'NF==9' | \
    awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_down.bed.temp - | \
    awk 'NF==10' | sort -k 2n | cut -f1-4 > reference_down.temp

awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_down.bed.temp reference_query_down_raw.bed | \
    awk 'NF==9' | \
    awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_down.bed.temp - | \
    awk 'NF==10' | sort -k 6n | cut -f5-8 > query_down.temp

# Create final downstream output
python2 "$script_dir/merge_bed_regions.py" reference_down.temp | paste sp1 - | cut -f1,3-4 > sp1.down.temp
python2 "$script_dir/merge_bed_regions.py" query_down.temp | paste sp2 - | cut -f1,3-4 | cat sp1.down.temp - > down.txt

rm sp1.down.temp

# Validate downstream coverage against length threshold
if [ -n "$length_threshold" ]; then
    awk -v a="$length3" -v b="$length4" -v len_threshold="$length_threshold" \
        'a == len_threshold && b == len_threshold' down.txt | \
        awk '{print $1"\t"$3-$2+1"\tdown"}' | \
        awk -va="$pair" -v len_threshold="$length_threshold" \
        '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_down
else
    # If no length threshold specified, create empty check file
    touch check_down
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downstream processing completed"

################################################################################
# SUMMARY
################################################################################
echo "=========================================="
echo "Processing Summary"
echo "=========================================="
echo "Pair: $pair"
echo "Distance threshold: $distance"
echo "Length threshold: ${length_threshold:-not specified}"
echo ""
echo "Output files generated:"
echo "  - up.txt (merged upstream regions)"
echo "  - down.txt (merged downstream regions)"
echo "  - check_up (upstream validation)"
echo "  - check_down (downstream validation)"
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] All processing completed successfully!"
echo "=========================================="