#!/bin/bash

# Check if all required parameters are provided
if [ $# -ne 5 ]; then
    echo "Usage: $0 <first_dir> <second_dir> <flank_size>"
    echo "Example: $0 5k 10k 10000 5 20"
    exit 1
fi

# Assign input parameters
first_dir=$1      # First directory name (e.g., 5k)
second_dir=$2     # Second directory name (e.g., 10k)
flank_size=$3     # Flanking sequence size (e.g., 10000)
cpu1=$4	# promer thread (e.g., 5)
cpu2=$5	# common processing script thread (e.g., 20)

# Process remaining gene pairs using the first directory results
cd $workspace/Dupgene/process/​​${first_dir}
rm -rf $workspace/Dupgene/process/​​${second_dir}
mkdir -p $workspace/Dupgene/process/​​${second_dir}

# Extract gene IDs for upstream and downstream analysis
awk '$3=="up"{print $4}' ${first_dir}_up_down_check | sort -u > $workspace/Dupgene/process/​​${second_dir}/${second_dir}_up_list
awk '$3=="down"{print $4}' ${first_dir}_up_down_check | sort -u > $workspace/Dupgene/process/​​${second_dir}/${second_dir}_down_list

cd $workspace/Dupgene/process/​​${second_dir}
cat ${second_dir}_up_list ${second_dir}_down_list | sort -u > ${second_dir}_list

## Extract flanking sequences (upstream/downstream)
# Extract upstream sequences
grep -wf ${second_dir}_list $workspace/DupGen_finder/Dupgene.pairs | awk '{print $2"\n"$3}' | sort -u | grep -wf - $workspace/Dupgene.bed | \
bedtools flank -i - -g $workspace/Dupgene/genome.fa.len -l $flank_size -r 0 -s > up.Dupgene.bed

# Extract downstream sequences
grep -wf ${second_dir}_list $workspace/DupGen_finder/Dupgene.pairs | awk '{print $2"\n"$3}' | sort -u | grep -wf - $workspace/Dupgene.bed | \
bedtools flank -i - -g $workspace/Dupgene/genome.fa.len -l 0 -r $flank_size -s > down.Dupgene.bed

## Extract FASTA sequences from genome
bedtools getfasta -s -fi $genome -bed up.Dupgene.bed -fo up.Dupgene.temp.fa -name
bedtools getfasta -s -fi $genome -bed down.Dupgene.bed -fo down.Dupgene.temp.fa -name

# Clean FASTA headers
sed 's/::/\t/' up.Dupgene.temp.fa | cut -f 1 > up.Dupgene.fa
sed 's/::/\t/' down.Dupgene.temp.fa | cut -f 1 > down.Dupgene.fa
rm up.Dupgene.temp.fa down.Dupgene.temp.fa

# Prepare temporary directory
rm -rf $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP
mkdir -p $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP
cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP

# Generate gene ID files
rm get_id.sh* 2>/dev/null
for i in $(cat ../${second_dir}_list)
do
echo "mkdir -p $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i};cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i};awk -va=\"$i\" '\$1==a {print \$2\"\\n\"\$3}' $workspace/DupGen_finder/Dupgene.pairs > ${i}.id"
done &> get_id.sh
ParaFly -c get_id.sh -CPU ${cpu2}

# Promer analysis for upstream sequences
rm promer.up.sh* 2>/dev/null
for i in $(cat ../${second_dir}_up_list)
do
echo "cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i}; \
sed -n '1p' ${i}.id > a.id; \
sed -n '2p' ${i}.id > b.id; \
seqkit grep -f a.id ../../up.Dupgene.fa > a_up_${second_dir}.fa; \
seqkit grep -f b.id ../../up.Dupgene.fa > b_up_${second_dir}.fa; \
promer --maxmatch a_up_${second_dir}.fa b_up_${second_dir}.fa -p up; \
delta-filter -i \$identity_threshold up.delta > up_filtered.delta; \
show-coords -Trcdl up_filtered.delta > up_filtered.coords; \
mummerplot --png up_filtered.delta -p up" 
done &> promer.up.sh
ParaFly -c promer.up.sh -CPU ${cpu1} > up_script.log 2>&1

# Promer analysis for downstream sequences
rm promer.down.sh* 2>/dev/null
for i in $(cat ../${second_dir}_down_list)
do
echo "cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i}; \
sed -n '1p' ${i}.id > a.id; \
sed -n '2p' ${i}.id > b.id; \
seqkit grep -f a.id ../../down.Dupgene.fa > a_down_${second_dir}.fa; \
seqkit grep -f b.id ../../down.Dupgene.fa > b_down_${second_dir}.fa; \
promer --maxmatch a_down_${second_dir}.fa b_down_${second_dir}.fa -p down; \
delta-filter -i \$identity_threshold down.delta > down_filtered.delta; \
show-coords -Trcdl down_filtered.delta > down_filtered.coords; \
mummerplot --png down_filtered.delta -p down"
done &> promer.down.sh
ParaFly -c promer.down.sh -CPU ${cpu1} > down_script.log 2>&1

# Process upstream results
cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP
rm run_up.sh* 2>/dev/null
for i in $(cut -f 1 ../${second_dir}_up_list)
do
echo "cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i};export workspace=$workspace;export script=$script;\$script/up.sh -p ${i} -d \$distance -l $flank_size 2>/dev/null"
done &> run_up.sh
ParaFly -c run_up.sh -CPU ${cpu2}

# Process downstream results
cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP
rm run_down.sh* 2>/dev/null
for i in $(cut -f 1 ../${second_dir}_down_list)
do
echo "cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP/${i};export workspace=$workspace;export script=$script;\$script/down.sh -p ${i} -d \$distance -l $flank_size 2>/dev/null"
done &> run_down.sh
ParaFly -c run_down.sh -CPU ${cpu2}

# Collect results
cd $workspace/Dupgene/process/​​${second_dir}/${second_dir}_TEMP
find . -name "check_up" | xargs cat | awk 'NF==4' > ../${second_dir}_up_check
find . -name "check_down" | xargs cat | awk 'NF==4' > ../${second_dir}_down_check        
cat ../${second_dir}_up_check ../${second_dir}_down_check > ../${second_dir}_up_down_check
