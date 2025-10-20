#!/bin/bash

#############################################################################
# Script: DupGen_analysis.sh
# Description: Identify and classify gene duplications using DupGen_finder
#
# This script performs the following steps:
# 1. Prepares GFF and protein files for DupGen_finder analysis
# 2. Runs DIAMOND BLAST to identify homologous genes
# 3. Classifies duplications into WGD, TD, PD, TRD, and DSD
# 4. Removes redundancy following priority: WGD > TD > PD > TRD > DSD
# 5. Generates summary statistics and output files
#############################################################################

# Function to display usage
usage() {
    cat << EOF
Usage: $0 -tp <target_pep> -tg <target_gff> -op <outgroup_pep> -og <outgroup_gff> -w <workspace>

Required Arguments:
    -tp    Path to target species protein sequences (FASTA format)
    -tg    Path to target species GFF annotation file
    -op    Path to outgroup species protein sequences (FASTA format)
    -og    Path to outgroup species GFF annotation file
    -w     Working directory for analysis output

Example:
    $0 -tp target.pep.fa -tg target.gff -op outgroup.pep.fa -og outgroup.gff -w ./workspace

Output:
    workspace/DupGen_finder/
        ├── data/                    # Input data files
        ├── raw/                     # Raw DupGen_finder results
        ├── DupGenes.sum            # Summary statistics
        ├── Dupgene.txt             # Duplicate gene list with types
        ├── Dupgene.bed             # Duplicate gene locations (BED format)
        ├── Dupgene.pairs           # Duplicate gene pairs
        └── *.genes                 # Gene lists by duplication type

Requirements:
    - DupGen_finder-unique.pl
    - DIAMOND
    - Standard Unix tools (awk, sed, cut, sort, etc.)

EOF
    exit 1
}

# Initialize variables
pep=""
gff=""
outgroup_pep=""
outgroup_gff=""
workspace=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -tp)
            pep="$2"
            shift 2
            ;;
        -tg)
            gff="$2"
            shift 2
            ;;
        -op)
            outgroup_pep="$2"
            shift 2
            ;;
        -og)
            outgroup_gff="$2"
            shift 2
            ;;
        -w)
            workspace="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$pep" ] || [ -z "$gff" ] || [ -z "$outgroup_pep" ] || [ -z "$outgroup_gff" ] || [ -z "$workspace" ]; then
    echo "Error: All required arguments must be provided"
    echo ""
    usage
fi

# Validate input files
if [ ! -f "$pep" ]; then
    echo "Error: Target protein file not found: $pep"
    exit 1
fi

if [ ! -f "$gff" ]; then
    echo "Error: Target GFF file not found: $gff"
    exit 1
fi

if [ ! -f "$outgroup_pep" ]; then
    echo "Error: Outgroup protein file not found: $outgroup_pep"
    exit 1
fi

if [ ! -f "$outgroup_gff" ]; then
    echo "Error: Outgroup GFF file not found: $outgroup_gff"
    exit 1
fi

# Create workspace directory
mkdir -p "$workspace/DupGen_finder/data"
cd "$workspace/DupGen_finder/data" || exit 1

echo "=========================================="
echo "Duplicate Gene Analysis Pipeline"
echo "=========================================="
echo "Start time: $(date)"
echo ""
echo "Parameters:"
echo "  Target protein:   $pep"
echo "  Target GFF:       $gff"
echo "  Outgroup protein: $outgroup_pep"
echo "  Outgroup GFF:     $outgroup_gff"
echo "  Workspace:        $workspace"
echo ""

#############################################################################
# Step 1: Prepare GFF files for DupGen_finder
#############################################################################
echo "[Step 1] Preparing GFF files..."

# Extract mRNA features from target GFF
cut -d ';' -f 1 "$gff" | awk '$3=="mRNA" {print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' > target.gff

# Extract mRNA features from outgroup GFF
cut -d ';' -f 1 "$outgroup_gff" | awk '$3=="mRNA" {print $1"\t"$9"\t"$4"\t"$5}' | sed 's/ID=//' > outgroup.gff

# Combine target and outgroup GFF files
cat target.gff outgroup.gff > target_outgroup.gff

echo "  - Target GFF: $(wc -l < target.gff) mRNA features"
echo "  - Outgroup GFF: $(wc -l < outgroup.gff) mRNA features"
echo ""

#############################################################################
# Step 2: Run DIAMOND BLAST for homology detection
#############################################################################
echo "[Step 2] Running DIAMOND BLAST..."

# Create DIAMOND database for target proteins
echo "  - Creating target protein database..."
diamond makedb --in "$pep" --db target.pep --quiet

# BLAST target proteins against themselves
echo "  - Running target self-BLAST..."
diamond blastp --query "$pep" --db target.pep --out target.blast \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 10 --threads 30 --quiet --sensitive

# Create DIAMOND database for outgroup proteins
echo "  - Creating outgroup protein database..."
diamond makedb --in "$outgroup_pep" --db outgroup.pep --quiet

# BLAST target proteins against outgroup proteins
echo "  - Running target vs outgroup BLAST..."
diamond blastp --query "$pep" --db outgroup.pep --out target_outgroup.blast \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 10 --threads 30 --quiet --sensitive

echo ""

#############################################################################
# Step 3: Run DupGen_finder to identify duplications
#############################################################################
echo "[Step 3] Running DupGen_finder..."

$DupGen_finder/DupGen_finder-unique.pl -i "$workspace/DupGen_finder/data" -t target -c outgroup -o ../raw -m 25 -s 5

echo "  - DupGen_finder analysis completed"
echo ""

#############################################################################
# Step 4: Remove redundancy in duplication categories
# Priority order: WGD > TD > PD > TRD > DSD
#############################################################################
echo "[Step 4] Removing redundancy from duplication categories..."

cd "$workspace/DupGen_finder/raw" || exit 1

# Extract WGD/Segmental Duplication genes
sed '1d' target.wgd.pairs-unique | cut -f 1,3 | sed 's/\t/\n/' | sort -u > ../wgd_or_sd.genes

# Extract genes from each duplication type to temporary files
sed '1d' target.tandem.pairs-unique | cut -f 1,3 | sed 's/\t/\n/' | sort -u > tandem.temp
sed '1d' target.proximal.pairs-unique | cut -f 1,3 | sed 's/\t/\n/' | sort -u > proximal.temp
sed '1d' target.transposed.pairs-unique | cut -f 1,3 | sed 's/\t/\n/' | sort -u > transposed.temp
sed '1d' target.dispersed.pairs-unique | cut -f 1,3 | sed 's/\t/\n/' | sort -u > dispersed.temp

# Remove tandem genes already classified as WGD/SD
cat ../wgd_or_sd.genes tandem.temp | sort | uniq -c | awk '$1==2{print $2}' > tandem.redundant
cat tandem.redundant tandem.temp | sort | uniq -c | awk '$1==1{print $2}' > ../tandem.genes

# Remove proximal genes already classified as WGD/SD or TD
cat ../wgd_or_sd.genes ../tandem.genes proximal.temp | sort | uniq -c | awk '$1==2{print $2}' > proximal.redundant
cat proximal.redundant proximal.temp | sort | uniq -c | awk '$1==1{print $2}' > ../proximal.genes

# Remove transposed genes already classified as WGD/SD, TD, or PD
cat ../wgd_or_sd.genes ../tandem.genes ../proximal.genes transposed.temp | sort | uniq -c | awk '$1==2{print $2}' > transposed.redundant
cat transposed.redundant transposed.temp | sort | uniq -c | awk '$1==1{print $2}' > ../transposed.genes

# Remove dispersed genes already classified in other categories
cat ../wgd_or_sd.genes ../tandem.genes ../proximal.genes ../transposed.genes dispersed.temp | sort | uniq -c | awk '$1==2{print $2}' > dispersed.redundant
cat dispersed.redundant dispersed.temp | sort | uniq -c | awk '$1==1{print $2}' > ../dispersed.genes

# Clean up temporary files
rm *.redundant *.temp

# Copy singleton genes
cp target.singletons ../singleton.genes

echo "  - Redundancy removal completed"
echo ""

#############################################################################
# Step 5: Generate summary statistics
#############################################################################
echo "[Step 5] Generating summary statistics..."

cd "$workspace/DupGen_finder" || exit 1

# Count genes in each category
wc -l wgd_or_sd.genes | awk '{print $2"\t"$1}' | sed 's/wgd_or_sd/WGD_or_SD/' > temp1
wc -l tandem.genes | awk '{print $2"\t"$1}' | sed 's/tandem/TD/' > temp2
wc -l proximal.genes | awk '{print $2"\t"$1}' | sed 's/proximal/PD/' > temp3
wc -l transposed.genes | awk '{print $2"\t"$1}' | sed 's/transposed/TRD/' > temp4
wc -l dispersed.genes | awk '{print $2"\t"$1}' | sed 's/dispersed/DSD/' > temp5
wc -l singleton.genes | awk '{print $2"\t"$1}' | sed 's/singleton/Singleton/' > temp6
wc -l ./data/target.gff | awk '{print "All\t"$1}' > temp7

# Calculate total number of duplicated genes
cat temp1 temp2 temp3 temp4 temp5 | awk -F "\t" '{sum += $2};END {print sum}' | awk '{print "Duplication\t"$1}' > temp8

# Create summary file
cat temp7 temp8 temp1 temp2 temp3 temp4 temp5 temp6 | sed 's/\.genes//' > DupGenes.sum

# Clean up temporary files
rm temp*

echo "  - Summary statistics saved to DupGenes.sum"
echo ""

#############################################################################
# Step 6: Create duplicate gene list with type annotations
#############################################################################
echo "[Step 6] Creating duplicate gene lists..."

cd "$workspace/DupGen_finder" || exit 1

# Combine all duplicate genes with their type labels
for i in tandem proximal transposed dispersed; do
    awk -va="$i" '{print a"\t"$0}' ${i}.genes
done | sed -e 's/tandem/TD/' -e 's/proximal/PD/' -e 's/transposed/TRD/' -e 's/dispersed/DSD/' > Dupgene.txt

echo "  - Duplicate gene list saved to Dupgene.txt"
echo ""

#############################################################################
# Step 7: Create BED file with duplicate gene locations
#############################################################################
echo "[Step 7] Creating BED file with gene locations..."

# Extract mRNA features and convert to BED format
# Note: BED format uses 0-based start coordinates, GFF uses 1-based
awk '$3=="mRNA"' "$gff" | cut -d ';' -f 1 | sed 's/ID=//' | \
    awk '{print $1"\t"($4-1)"\t"$5"\t"$9"\t"$7}' | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$4]}' Dupgene.txt - | \
    awk 'NF==6{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' > Dupgene.bed

echo "  - BED file saved to Dupgene.bed"
echo ""

#############################################################################
# Step 8: Create duplicate gene pairs file
#############################################################################
echo "[Step 8] Creating duplicate gene pairs..."

cd "$workspace/DupGen_finder" || exit 1

# Process tandem duplications
awk '$1=="TD"' Dupgene.txt > Temp_TD
sed '1d' ./raw/target.tandem.pairs-unique | cut -f 1,3 | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$1]}' Temp_TD - | awk 'NF==3' | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$2]}' Temp_TD - | \
    awk 'NF==4 {print "TD"NR"\t"$1"\t"$2"\tTD"}' > Temp_1

# Process proximal duplications
awk '$1=="PD"' Dupgene.txt > Temp_PD
sed '1d' ./raw/target.proximal.pairs-unique | cut -f 1,3 | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$1]}' Temp_PD - | awk 'NF==3' | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$2]}' Temp_PD - | \
    awk 'NF==4 {print "PD"NR"\t"$1"\t"$2"\tPD"}' > Temp_2

# Process transposed duplications
awk '$1=="TRD"' Dupgene.txt > Temp_TRD
sed '1d' ./raw/target.transposed.pairs-unique | cut -f 1,3 | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$1]}' Temp_TRD - | awk 'NF==3' | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$2]}' Temp_TRD - | \
    awk 'NF==4 {print "TRD"NR"\t"$1"\t"$2"\tTRD"}' > Temp_3

# Process dispersed duplications
awk '$1=="DSD"' Dupgene.txt > Temp_DSD
sed '1d' ./raw/target.dispersed.pairs-unique | cut -f 1,3 | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$1]}' Temp_DSD - | awk 'NF==3' | \
    awk 'NR==FNR{a[$2]=$1;next}{print $0"\t"a[$2]}' Temp_DSD - | \
    awk 'NF==4 {print "DSD"NR"\t"$1"\t"$2"\tDSD"}' > Temp_4

# Combine all pairs
cat Temp_1 Temp_2 Temp_3 Temp_4 > Dupgene.pairs

# Clean up temporary files
rm Temp_*

echo "  - Duplicate gene pairs saved to Dupgene.pairs"
echo ""

#############################################################################
# Final summary
#############################################################################
echo "=========================================="
echo "Analysis completed successfully!"
echo "=========================================="
echo "End time: $(date)"
echo ""
echo "Output directory: $workspace/DupGen_finder/"
echo ""
echo "Summary of results:"
cat "$workspace/DupGen_finder/DupGenes.sum"
echo ""
echo "Main output files:"
echo "  - DupGenes.sum     : Summary statistics"
echo "  - Dupgene.txt      : Duplicate genes with type labels"
echo "  - Dupgene.bed      : Duplicate gene locations (BED format)"
echo "  - Dupgene.pairs    : Duplicate gene pairs"
echo "  - *.genes          : Individual gene lists by duplication type"
echo ""
echo "Note: Gene count and gene pair count may not correspond 1:1"
echo "      due to one-to-many relationships in duplications."
echo ""