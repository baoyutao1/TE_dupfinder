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

###### Upstream Processing ######
# Generate upstream temporary files
sed '1,4d' up_filtered.coords |awk '{print "chr\t"$1"\t"$2"\t"$16}' | awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 1
sed '1,4d' up_filtered.coords |awk '{print "chr\t"$3"\t"$4"\t"$17}' | awk '$2<$3{print $0}$2>$3{print $1"\t"$3"\t"$2"\t"$4}' > 2
paste 1 2|sort -k 7nr -k 3nr > reference_query_up_raw.bed
paste 1 2|sort -k 3nr -k 7nr > reference_query_up_raw_sort.bed
cut -f 4 1 | sort -u > sp1
cut -f 4 2 | sort -u > sp2

# Get gene lengths from Dupgene.bed
length1=$(grep -wf sp1 ../../up.Dupgene.bed | awk '{print $3-$2}')
length2=$(grep -wf sp2 ../../up.Dupgene.bed | awk '{print $3-$2}')

# Keep rows where the normalized number of column 2 to column 6 values (relative to previous row) 
python2 $script/filter_raw_data.py reference_query_up_raw.bed -d -500 -u 500 > temp1
python2 $script/filter_raw_data.py reference_query_up_raw_sort.bed -d -500 -u 500 > temp2
cat temp1 temp2|sort -u|sort -k 7nr -k 3nr > reference_query_up_clean.bed
#awk -va=$length1 -vb=$length2 '{ratio = (b-$7+1)/(a-$3+1);if (ratio >= 0.5 && ratio <= 2) print $0}' reference_query_up_raw.bed > reference_query_up_clean.bed
cut -f 1-4 reference_query_up_clean.bed | sort -k 3nr -k 2nr > reference_up_raw.bed
cut -f 5-8 reference_query_up_clean.bed | sort -k 3nr -k 2nr > query_up_raw.bed
rm 1 2 temp1 temp2 reference_query_up_raw_sort.bed

# Filter BED files: remove lines completely contained within previous lines
python2 $script/filter_bed.py reference_up_raw.bed > reference_up_clean.bed
python2 $script/filter_bed.py query_up_raw.bed > query_up_clean.bed

# Process upstream regions with distance threshold
python2 $script/up.py -c $length1 -d $distance reference_up_clean.bed > reference_up.bed.temp
python2 $script/up.py -c $length2 -d $distance query_up_clean.bed > query_up.bed.temp

# Merge and filter results
awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_up.bed.temp reference_query_up_raw.bed |\
awk 'NF==9' | awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_up.bed.temp - |\
awk 'NF==10' | sort -k 3nr | cut -f1-4 > reference_up.temp

awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$1$2$3$4]}' reference_up.bed.temp reference_query_up_raw.bed |\
awk 'NF==9' | awk 'NR==FNR{a[$1$2$3$4]=$1;next}{print $0"\t"a[$5$6$7$8]}' query_up.bed.temp - |\
awk 'NF==10' | sort -k 7nr | cut -f5-8 > query_up.temp

# Merge regions and create final upstream file
python2 $script/merge_bed_regions.py reference_up.temp | paste sp1 - | cut -f1,3-4 > sp1.up.temp
python2 $script/merge_bed_regions.py query_up.temp | paste sp2 - | cut -f1,3-4 | cat sp1.up.temp - > up.txt
rm sp1.up.temp

# Check if upstream coverage meets threshold (90% of custom length_threshold)
#awk '{print $1"\t"$3-$2+1"\tup"}' up.txt | awk -va="$pair" -v len_threshold="$length_threshold" '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_up
awk -v a="$length1" -v b="$length2" -v len_threshold="$length_threshold" 'a == len_threshold && b == len_threshold' up.txt|\
awk '{print $1"\t"$3-$2+1"\tup"}' | awk -va="$pair" -v len_threshold="$length_threshold" '$2 > len_threshold * 0.9 {print $0"\t"a}' > check_up

echo "Upstream processing completed for pair: $pair with distance: $distance and length threshold: $length_threshold"