#!/bin/bash
f1=$1 # psicem bed file
f2=$2 # chromap bed file
out=$3 # output directory
cut -f 4 $f1 | sort | uniq | sort > barcode_piscem
cut -f 4 $f2 | sort | uniq | sort > barcode_chromap
comm -13 barcode_piscem barcode_chromap > ${out}/missing_barcode_piscem.txt
comm -23 barcode_piscem barcode_chromap > ${out}/missing_barcode_chromap.txt
rm barcode_piscem barcode_chromap