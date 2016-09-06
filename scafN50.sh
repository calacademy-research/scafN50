#!/bin/bash
# run All, 201, 300, 500, and 10,000, N50 analysis for the scaffold/contigs in the scafSeq file of arg 1
# results are output to scafN50s.txt
# 26Jan15 add 300
if [ $# -ne 1 ]; then
	echo
    echo "    scafN50.sh computes the scaffold N50 for the scafSeq file of argument 1."
	echo "    It computes N50s for all the scaffolds/contigs and also excluding scaffold of several sizes and below."
	echo "    The output is in a file named scafN50s.txt"
	echo
else
	out=scafN50s.txt
	scaflens.py $1 -n -C 201 300 500 1000 10000 > $out
fi
