#!/bin/bash

# Check the number of arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <variable>"
    exit 1
fi

variable=$1

# 处理down数据
{
    grep -wf sp1 down.txt |  awk 'NF==3' | awk -vb=${variable} '{print $1"\t"$2"\t"$3"\tdown\t"b}' 
    grep -wf sp2 down.txt |  awk 'NF==3' | awk -vb=${variable} '{print $1"\t"$2"\t"$3"\tdown\t"b}' 
} > down_tract.txt

echo "Processing completed. Results saved in $variable:down_tract.txt"