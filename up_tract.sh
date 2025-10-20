#!/bin/bash

# Check the number of arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <variable>"
    exit 1
fi

variable=$1

# Calculate gene lengths from BED file
length1=$(grep -wf sp1 "../../up.Dupgene.bed" | awk '{print $3-$2}')
length2=$(grep -wf sp2 "../../up.Dupgene.bed" | awk '{print $3-$2}')

# Process sp1 data: 
# 1. Find matching entries in up.txt
# 2. Filter lines with exactly 3 fields
# 3. Calculate reverse positions using gene length
# 4. Add annotation type and variable value
grep -wf sp1 up.txt | awk 'NF==3' | awk -v len="$length1" -v var="$variable" '{print $1"\t"len-$2"\t"len-$3"\tup\t"var}' > temp

# Process sp2 data and combine with sp1 results
grep -wf sp2 up.txt | awk 'NF==3' | awk -v len="$length2" -v var="$variable" '{print $1"\t"len-$2"\t"len-$3"\tup\t"var}' | cat temp - > up_tract.txt

# Rename temporary file to final result
rm temp

echo "Processing completed. Results saved in $variable:up_tract.txt"