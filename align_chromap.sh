#!/usr/bin/env bash

#SBATCH --partition=cbcb
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH --qos=high
#SBATCH --account=cbcb
#SBATCH --constraint=EPYC-7313

ulimit -n 2048
ref_file=/fs/cbcb-lab/rob/students/noor/Uncertainity/brain_sim/annotation/hg38.fa.gz
index_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/chromap_analysis/ind_hg38
#out_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/chromap_analysis
out_dir=/fs/nexus-projects/scATAC-seq/chromap

out_suff=$1
barcode_file=$2
read1_file=$3
read2_file=$4
whitelist_file=/fs/cbcb-lab/rob/students/noor/Atacseq/cellranger-atac-2.1.0/lib/python/atac/barcodes/737K-cratac-v1.txt

if [[ -n $4 ]]
    then
        /usr/bin/time -o ${out_dir}/time_align_${out_suff}.out /fs/cbcb-lab/rob/students/noor/chromap2/chromap/chromap \
            --preset atac \
            -t 16 \
            -x $index_dir \
            -r $ref_file \
            -1 $read1_file \
            -2 $read2_file \
            -b $barcode_file \
	    -o ${out_dir}/map_${out_stuff}.bed \
	    --summary A_summary \
	    --barcode-whitelist ${whitelist_file} \
	    --read-format bc:0:-1:-
            #-o ${out_dir}/map_${out_suff}.bed \
#	    --SAM \
#            -o ${out_dir}/map_${out_suff}.SAM \
else
        /usr/bin/time -o ${out_dir}/time_align_${out_suff}.out /fs/cbcb-lab/rob/students/noor/chromap2/chromap \
            --preset atac \
            -t 16 \
            -x $index_dir \
            -r $ref_file \
            -1 $read1_file \
            -b $barcode_file \
	    --SAM \
            -o ${out_dir}/map_out_${out_suff}.SAM \
	    --summary A_summary
fi

## Example
## sbatch helper_scripts/align_chromap.sh samp1285 sub_seq/barcode_subsamp.fq sub_seq/read1_subsamp.fq sub_seq/read2_subsamp.fq
