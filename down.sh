#!/bin/bash

# Initialize variables with default values
distance=500
pair=""
length_threshold=""

# Parse command-line options
while getopts "d:p:l:" opt; do
    case $opt in
        d) distance="$OPTARG" ;;
        p) pair="$OPTARG" ;;
        l) length_threshold="$OPTARG" ;;  # New option for length threshold
        *) echo "Usage: $0 -p <pair> [-d <distance>] [-l <length_threshold>]"; exit 1 ;;
    esac
done

# Check if pair parameter exists
if [ -z "$pair" ]; then
    echo "Usage: $0 -p <pair> [-d <distance>] [-l <length_threshold>]"
    exit 1
fi

# Create directory and enter
base_dir=$(pwd) 
cd "$base_dir" || exit 1

###### Downstream Processing ######
# Generate downstream temporary files
sed '1,4d' down_filtered.coords |awk '{print "chr\t"$1"\t"$2"\t"$16}' | awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 1
sed '1,4d' down_filtered.coords |awk '{print "chr\t"$3"\t"$4"\t"$17}' | awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 2
paste 1 2 |sort -k 6n -k 2n > reference_query_down_raw.bed
paste 1 2 |sort -k 2n -k 6n > reference_query_down_raw_sort.bed
cut -f 4 1 | sort -u > sp1
cut -f 4 2 | sort -u > sp2

# Get gene lengths from Dupgene.bed
length3=$(grep -wf sp1 ../../down.Dupgene.bed | awk '{print $3-$2}')
length4=$(grep -wf sp2 ../../down.Dupgene.bed | awk '{print $3-$2}')

# Keep rows where the normalized number of column 2 to column 6 values (relative to previous row) 
python2 $script/filter_raw_data.py reference_query_down_raw.bed -d -500 -u 500 > temp1
python2 $script/filter_raw_data.py reference_query_down_raw_sort.bed -d -500 -u 500 > temp2
cat temp1 temp2|sort -u|sort -k 6n -k 2n > reference_query_down_clean.bed
#awk '$2 != 0 {ratio = ($6+1)/($2+1);if (ratio >= 0.5 && ratio <= 2) print $0}' reference_query_down_raw.bed > reference_query_down_clean.bed
cut -f 1-4 reference_query_down_clean.bed | sort -k 2n -k 3n > reference_down_raw.bed
cut -f 5-8 reference_query_down_clean.bed | sort -k 2n -k 3n > query_down_raw.bed
rm 1 2 temp1 temp2 reference_query_down_raw_sort.bed

# Filter BED files for downstream
python2 $script/filter_bed.py reference_down_raw.bed > reference_down_clean.bed
python2 $script/filter_bed.py query_down_raw.bed > query_down_clean.bed

# Process downstream regions
python2 $script/down.py reference_down_clean.bed -d $distance > reference_down.bed.temp
python2 $script/down.py query_down_clean.bed -d $distance > query_down.bed.temp

# Merge and filter downstream results
awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_down.bed.temp reference_query_down_raw.bed |\
awk 'NF==9' | awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_down.bed.temp - |\
awk 'NF==10' | sort -k 2n | cut -f1-4 > reference_down.temp

awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_down.bed.temp reference_query_down_raw.bed |\
awk 'NF==9' | awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_down.bed.temp - |\
awk 'NF==10' | sort -k 6n | cut -f5-8 > query_down.temp

# Merge regions and create final downstream file
python2 $script/merge_bed_regions.py reference_down.temp | paste sp1 - | cut -f1,3-4 > sp1.down.temp
python2 $script/merge_bed_regions.py query_down.temp | paste sp2 - | cut -f1,3-4 | cat sp1.down.temp - > down.txt
rm sp1.down.temp

# Check if downstream coverage meets threshold (90% of custom length_threshold)
#awk '{print $1"\t"$3-$2+1"\tdown"}' down.txt | awk -va="$pair" -v len_threshold="$length_threshold" '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_down
awk -v a="$length3" -v b="$length4" -v len_threshold="$length_threshold" 'a == len_threshold && b == len_threshold' down.txt|\
awk '{print $1"\t"$3-$2+1"\tdown"}' | awk -va="$pair" -v len_threshold="$length_threshold" '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_down

echo "Downstream processing completed for pair: $pair with distance: $distance and length threshold: $length_threshold"