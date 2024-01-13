#!/usr/bin/env bash

#SBATCH --partition=cbcb
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH --qos=high
#SBATCH --account=cbcb

ulimit -n 2048
ref_file=/fs/cbcb-lab/rob/students/noor/Uncertainity/brain_sim/annotation/hg38.fa.gz
index_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/chromap_analysis/ind_hg38
out_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/chromap_analysis
out_suff=$1
barcode_file=$2
read1_file=$3
read2_file=$4

if [[ -n $4 ]]
    then
        /usr/bin/time -o ${out_dir}/time_align_${out_suff}.out /fs/cbcb-lab/rob/students/noor/chromap/chromap \
            --preset atac \
            -t 16 \
            -x $index_dir \
            -r $ref_file \
            -1 $read1_file \
            -2 $read2_file \
            -b $barcode_file \
            -o ${out_dir}/map_out_${out_suff}.bed
else
        /usr/bin/time -o ${out_dir}/time_align_${out_suff}.out /fs/cbcb-lab/rob/students/noor/chromap/chromap \
            --preset atac \
            -t 16 \
            -x $index_dir \
            -r $ref_file \
            -1 $read1_file \
            -b $barcode_file \
            -o ${out_dir}/map_out_${out_suff}.bed
fi