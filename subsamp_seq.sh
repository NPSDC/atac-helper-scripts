#!/bin/bash

fastq_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fastqs
base_filename=10k_pbmc_ATACv2_nextgem_Chromium_Controller_S3_L001
seqkit_path=/fs/cbcb-lab/rob/students/noor/Atacseq/seqkit 
out_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/$1
seed=$2

read1_file=${out_dir}/read1_subsamp.fq
read2_file=${out_dir}/read2_subsamp.fq
barcode_file=${out_dir}/barcode_subsamp.fq

read1_temp=${out_dir}/read1_temp.fq
read2_temp=${out_dir}/read2_temp.fq
barcode_temp=${out_dir}/barcode_temp.fq

rem_dup_bar=true

$seqkit_path sample -p 0.000005 -s $seed $fastq_dir/${base_filename}_R1_001.fastq.gz -j 16  -o $read1_file
$seqkit_path sample -p 0.000005 -s $seed $fastq_dir/${base_filename}_R3_001.fastq.gz -j 16  -o $read2_file
$seqkit_path sample -p 0.000005 -s $seed $fastq_dir/${base_filename}_R2_001.fastq.gz -j 16  -o $barcode_file

### We can have barcodes repeated ideally (each barcode corresponds to cell)

### Remove barcodes that appear more than once (just keep the first entry) and accordingly update
### the fastq files
if [ $rem_dup_bar ]; then
    > ${read1_temp}
    > ${read2_temp}
    > ${barcode_temp}
    mapfile -t barcodes < <( sed -n '2~4p' ${barcode_file} | sort | uniq)
    for barcode in ${barcodes[@]};do
        mapfile -t header_barcodes < <(grep -B 1 $barcode ${barcode_file} | grep -v -- -- | awk 'NR %2 != 0')
        header_barcode=${header_barcodes[0]}
        header_f1=$(echo ${header_barcode} | sed 's/ 2:N:0/ 1:N:0/')
        grep -h -A 3 "${header_f1}" ${read1_file} >> ${read1_temp}
        grep -h -A 3 "${header_barcode}" ${barcode_file} >> ${barcode_temp}
        header_f2=$(echo ${header_barcode} | sed 's/ 2:N:0/ 3:N:0/')
        grep -h -A 3 "${header_f2}" ${read2_file} >> ${read2_temp}
    done
    mv $read1_temp $read1_file
    mv $read2_temp $read2_file
    mv $barcode_temp $barcode_file
fi
