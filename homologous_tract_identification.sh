#!/bin/bash

##############################################################################
# Script: homologous_tract_identification.sh
# Description: Analyze duplicated genes and identify homologous regions in 
#              upstream and downstream flanking sequences through iterative 
#              extension rounds
##############################################################################

set -e  # Exit on error

# Default values
workspace=""
extended_length=""
genome=""
general_thread=20
identity_threshold=60
promer_thread=5
distance=500

# Help message
show_help() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Description:
    This script analyzes duplicated genes by extracting and comparing upstream 
    and downstream flanking sequences. It performs iterative rounds of homology 
    searches using MUMmer/promer to identify conserved syntenic regions.

Required Options:
    -w, --workspace DIR         Working directory path
    -el, --extended-length NUM  Extended length parameter for sequence generation
    -g, --genome FILE           Reference genome FASTA file path

Optional Parameters:
    -gt, --general-thread NUM   Number of threads for general parallel tasks (default: 20)
    -it, --identity-threshold NUM  Identity threshold for delta-filter (default: 60)
    -pt, --promer-thread NUM    Number of threads for promer tasks (default: 5)
    -d, --distance NUM          Distance parameter for analysis (default: 500)
    -h, --help                  Display this help message and exit

Required Input Files:
    - \$workspace/DupGen_finder/Dupgene.bed
    - \$workspace/DupGen_finder/Dupgene.pairs
    - Helper scripts in \$script directory

Output Files:
    - \$workspace/Dupgene/Dupgene.up.bed
    - \$workspace/Dupgene/Dupgene.down.bed
    - \$workspace/Dupgene/Dupgene.up_down.bed
    - \$workspace/Dupgene/Dupgene.up_down.breakpoint.txt

Dependencies:
    - samtools, bedtools, seqkit, promer, delta-filter, show-coords, mummerplot, ParaFly

Example:
    $(basename "$0") -w /path/to/workspace -g /path/to/genome.fa -el 1000000 -gt 20 -it 70 -pt 5 -d 500

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -w|--workspace)
            workspace="$2"
            shift 2
            ;;
        -el|--extended-length)
            extended_length="$2"
            shift 2
            ;;
        -g|--genome)
            genome="$2"
            shift 2
            ;;
        -gt|--general-thread)
            general_thread="$2"
            shift 2
            ;;
        -it|--identity-threshold)
            identity_threshold="$2"
            shift 2
            ;;
        -pt|--promer-thread)
            promer_thread="$2"
            shift 2
            ;;
        -d|--distance)
            distance="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        *)
            echo "Error: Unknown option $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Validate required parameters
if [[ -z "$workspace" ]]; then
    echo "Error: Workspace directory (-w) is required"
    exit 1
fi

if [[ -z "$extended_length" ]]; then
    echo "Error: Extended length (-el) is required"
    exit 1
fi

if [[ -z "$genome" ]]; then
    echo "Error: Genome file (-g) is required"
    exit 1
fi

# Check if genome file exists
if [[ ! -f "$genome" ]]; then
    echo "Error: Genome file does not exist: $genome"
    exit 1
fi

# Check if required environment variable is set
if [[ -z "$script" ]]; then
    echo "Error: \$script environment variable must be set to scripts directory"
    exit 1
fi

# Create directory structure
echo "Creating directory structure..."
mkdir -p "$workspace/Dupgene/process"

# Index genome and create length file
echo "Indexing genome and creating length file..."
samtools faidx "$genome" -o "$workspace/Dupgene/genome.fa.fai"
cut -f 1,2 "$workspace/Dupgene/genome.fa.fai" > "$workspace/Dupgene/genome.fa.len"

# Generate number sequences
cd "$workspace/Dupgene/process"
echo "Generating extended length sequences..."
"$script/number_sequence_generator.py" -n "$extended_length" > number_sequence.txt

# Read sequences into round variables
for i in {1..10}; do
    declare "round$i=$(sed -n "${i}p" number_sequence.txt)"
done

# Process Round 1
echo "Processing Round 1..."
rm -rf "$workspace/Dupgene/process/round1"
mkdir -p "$workspace/Dupgene/process/round1"
cd "$workspace/Dupgene/process/round1"

# Extract upstream and downstream regions
echo "Extracting upstream regions (length: $round1)..."
bedtools flank -i "$workspace/DupGen_finder/Dupgene.bed" -g "$workspace/Dupgene/genome.fa.len" -l "$round1" -r 0 -s > up.Dupgene.bed

echo "Extracting downstream regions (length: $round1)..."
bedtools flank -i "$workspace/DupGen_finder/Dupgene.bed" -g "$workspace/Dupgene/genome.fa.len" -l 0 -r "$round1" -s > down.Dupgene.bed

# Extract sequences from genome
echo "Extracting sequences from genome..."
bedtools getfasta -s -fi "$genome" -bed up.Dupgene.bed -fo up.Dupgene.temp.fa -name
bedtools getfasta -s -fi "$genome" -bed down.Dupgene.bed -fo down.Dupgene.temp.fa -name

# Clean up sequence headers
sed 's/::/\t/' up.Dupgene.temp.fa | cut -f 1 > up.Dupgene.fa
sed 's/::/\t/' down.Dupgene.temp.fa | cut -f 1 > down.Dupgene.fa
rm up.Dupgene.temp.fa down.Dupgene.temp.fa

# Extract gene pair list
cut -f 1 "$workspace/DupGen_finder/Dupgene.pairs" > raw_list

# Create temporary working directory
rm -rf "$workspace/Dupgene/process/round1/round1_TEMP"
mkdir -p "$workspace/Dupgene/process/round1/round1_TEMP"
cd "$workspace/Dupgene/process/round1/round1_TEMP"

# Generate script to extract gene IDs
echo "Generating gene ID extraction script..."
rm -f get_id.sh
while IFS= read -r i; do
    echo "mkdir -p $workspace/Dupgene/process/round1/round1_TEMP/${i}; cd $workspace/Dupgene/process/round1/round1_TEMP/${i}; awk -va=\"$i\" '\$1==a {print \$2\"\\n\"\$3}' $workspace/DupGen_finder/Dupgene.pairs > ${i}.id"
done < ../raw_list > get_id.sh
ParaFly -c get_id.sh -CPU "$general_thread"

# Generate promer comparison script with better error handling
echo "Generating promer alignment script..."
rm -f promer.sh

while IFS= read -r i; do
    cat >> promer.sh << PROMER_CMD
cd $workspace/Dupgene/process/round1/round1_TEMP/${i}; \
sed -n '1p' ${i}.id > a.id; \
sed -n '2p' ${i}.id > b.id; \
seqkit grep -f a.id ../../up.Dupgene.fa > a_up_round1.fa; \
seqkit grep -f b.id ../../up.Dupgene.fa > b_up_round1.fa; \
\$MUMmer/promer --maxmatch a_up_round1.fa b_up_round1.fa -p up; \
\$MUMmer/delta-filter -i $identity_threshold up.delta > up_filtered.delta; \
\$MUMmer/show-coords -Trcdl up_filtered.delta > up_filtered.coords; \
\$MUMmer/mummerplot --png up_filtered.delta -p up; \
seqkit grep -f a.id ../../down.Dupgene.fa > a_down_round1.fa; \
seqkit grep -f b.id ../../down.Dupgene.fa > b_down_round1.fa; \
\$MUMmer/promer --maxmatch a_down_round1.fa b_down_round1.fa -p down; \
\$MUMmer/delta-filter -i $identity_threshold down.delta > down_filtered.delta; \
\$MUMmer/show-coords -Trcdl down_filtered.delta > down_filtered.coords; \
\$MUMmer/mummerplot --png down_filtered.delta -p down
PROMER_CMD
done < ../raw_list

echo "Running promer alignments (this may take a while)..."
echo "Total commands: $(wc -l < promer.sh)"
ParaFly -c promer.sh -CPU "$promer_thread" > script.log 2>&1 || true

# Identify directories with successful alignments (look for coords files instead of png)
echo "Collecting successful alignment results..."
cd "$workspace/Dupgene/process/round1/round1_TEMP"
find . -type f -name "*.gp" -exec dirname {} \; | sort -u | awk -F/ '{print $2}' | uniq > ../round1_list

echo "Successful alignments: $(wc -l < ../round1_list)"

# Run round1 analysis script
echo "Running round1 analysis..."
rm -f run.sh
while IFS= read -r i; do
    echo "cd $workspace/Dupgene/process/round1/round1_TEMP/${i}; export workspace=$workspace; export script=$script; \$script/round1.sh -p ${i} -d $distance -l $round1 2>/dev/null"
done < ../round1_list > run.sh
ParaFly -c run.sh -CPU "$general_thread"

# Collect check results
cd "$workspace/Dupgene/process/round1/round1_TEMP"
find . -name "check_up" | xargs cat > ../round1_check_up 2>/dev/null || touch ../round1_check_up
find . -name "check_down" | xargs cat | cat ../round1_check_up - | awk 'NF==4' > ../round1_up_down_check
rm -f ../round1_check_up

# Iterative extension rounds (2-10)
echo "Running iterative extension rounds 2-10..."
cd "$workspace/Dupgene/process"
export workspace=$workspace
"$script/extend_survey.sh" -fd round1 -sd round2 -f "$round2" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round2 -sd round3 -f "$round3" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round3 -sd round4 -f "$round4" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round4 -sd round5 -f "$round5" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round5 -sd round6 -f "$round6" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round6 -sd round7 -f "$round7" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round7 -sd round8 -f "$round8" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round8 -sd round9 -f "$round9" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"
"$script/extend_survey.sh" -fd round9 -sd round10 -f "$round10" -it "$identity_threshold" -pt "$promer_thread" -gt "$general_thread" -d "$distance"

# Extract filter IDs
echo "Extracting filter IDs..."
cd "$workspace/Dupgene/process"
cut -f 4 "$workspace/Dupgene/process/round10/round10_up_down_check" | sort -u > filter.id

# Summarize information across all rounds
echo "Summarizing information across all rounds..."
cd "$workspace/Dupgene/process"

# Generate concatenation scripts for all rounds
for j in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    rm -f "${j}_cat.sh"
    while IFS= read -r i; do
        cat >> "${j}_cat.sh" << CAT_CMD
cd $workspace/Dupgene/process/${j}/${j}_TEMP/${i}; \
awk -va=${i} '{print \$0"\tup\t"a}' up.txt > up.temp; \
awk -va=${i} '{print \$0"\tdown\t"a}' down.txt | cat up.temp - > up_down; \
rm -f up.temp
CAT_CMD
    done < "$workspace/Dupgene/process/${j}/${j}_list"
done

# Execute concatenation scripts
for i in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    ParaFly -c "${i}_cat.sh" -CPU "$general_thread" > script.log 2>&1
done
rm -f FailedCommands script.log *.completed *_cat.sh

# Generate tract lists for each round
echo "Generating tract lists..."
for i in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    find "./${i}/${i}_TEMP" -name "up_down" | xargs cat | awk 'NF==5 && $4=="up" {print $5}' | sort -u > "${i}_up_list"
    find "./${i}/${i}_TEMP" -name "up_down" | xargs cat | awk 'NF==5 && $4=="down" {print $5}' | sort -u > "${i}_down_list"
    awk '$3=="up"{print $4}' "./${i}/${i}_up_down_check" | sort -u | grep -wvf - "${i}_up_list" | sort -u > "${i}_up.tract.list"
    awk '$3=="down"{print $4}' "./${i}/${i}_up_down_check" | sort -u | grep -wvf - "${i}_down_list" | sort -u > "${i}_down.tract.list"
    rm -f "${i}_up_list" "${i}_down_list"
done

# Create summary tract lists
for i in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    awk -va="${i}" '{print a"\t"$0}' "${i}_up.tract.list"
done > up.tract.list.sum

for i in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    awk -va="${i}" '{print a"\t"$0}' "${i}_down.tract.list"
done > down.tract.list.sum

# Process upstream tracts
echo "Processing upstream tracts..."
for j in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    rm -f "${j}_up.tract.sh"
    while IFS= read -r i; do
        echo "cd $workspace/Dupgene/process/${j}/${j}_TEMP/${i}; $script/up_tract.sh ${i}"
    done < "${j}_up.tract.list" > "${j}_up.tract.sh"
done

for j in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    ParaFly -c "${j}_up.tract.sh" -CPU "$general_thread"
done

find . -name 'up_tract.txt' | xargs cat | grep -wvf filter.id - > ../Dupgene.up.txt

# Process downstream tracts
echo "Processing downstream tracts..."
for j in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    rm -f "${j}_down.tract.sh"
    while IFS= read -r i; do
        echo "cd $workspace/Dupgene/process/${j}/${j}_TEMP/${i}; $script/down_tract.sh ${i}"
    done < "${j}_down.tract.list" > "${j}_down.tract.sh"
done

for j in round1 round2 round3 round4 round5 round6 round7 round8 round9 round10; do
    ParaFly -c "${j}_down.tract.sh" -CPU "$general_thread"
done

find . -name 'down_tract.txt' | xargs cat | grep -wvf filter.id - > ../Dupgene.down.txt

rm -f *tract.list

# Generate final BED files and breakpoint coordinates
echo "Generating final output files..."
cd "$workspace/Dupgene"

# Process upstream regions
awk 'NR==FNR{a[$4]=$1"\t"$2"\t"$3"\t"$6;next}{print $0"\t"a[$1]}' \
    "$workspace/DupGen_finder/Dupgene.bed" Dupgene.up.txt | \
    awk '$9=="+"{print $6"\t"$7-$2"\t"$7-$3"\t"$1"\t"$4"\t"$9"\t"$5}
         $9=="-"{print $6"\t"$8+$3"\t"$8+$2"\t"$1"\t"$4"\t"$9"\t"$5}' > Dupgene.up.bed

awk '$6=="+"{print $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}
     $6=="-"{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' Dupgene.up.bed > Dupgene.up.breakpoint.txt

# Process downstream regions
awk 'NR==FNR{a[$4]=$1"\t"$2"\t"$3"\t"$6;next}{print $0"\t"a[$1]}' \
    "$workspace/DupGen_finder/Dupgene.bed" Dupgene.down.txt | \
    awk '$9=="+"{print $6"\t"$8+$2"\t"$8+$3"\t"$1"\t"$4"\t"$9"\t"$5}
         $9=="-"{print $6"\t"$7-$2"\t"$7-$3"\t"$1"\t"$4"\t"$9"\t"$5}' > Dupgene.down.bed

awk '$6=="+"{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}
     $6=="-"{print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' Dupgene.down.bed > Dupgene.down.breakpoint.txt

# Combine and sort results
cat Dupgene.up.bed Dupgene.down.bed | sort -k 7,7V -k 5,5r > Dupgene.up_down.bed
cat Dupgene.up.breakpoint.txt Dupgene.down.breakpoint.txt | sort -k 7,7V -k 5,5r > Dupgene.up_down.breakpoint.txt

echo "============================================"
echo "Analysis completed successfully!"
echo "Output files:"
echo "  - $workspace/Dupgene/Dupgene.up.bed"
echo "  - $workspace/Dupgene/Dupgene.down.bed"
echo "  - $workspace/Dupgene/Dupgene.up_down.bed"
echo "  - $workspace/Dupgene/Dupgene.up_down.breakpoint.txt"
echo "============================================"