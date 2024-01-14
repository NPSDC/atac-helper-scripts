#!/usr/bin/env bash

#SBATCH --partition=cbcb
#SBATCH --cpus-per-task=16
#SBATCH --time=05:00:00
#SBATCH --qos=high
#SBATCH --account=cbcb

module load gcc/11.2.0
module load cmake/3.22.1
module load jemalloc/5.2.1

set -euo pipefail

ulimit -n 2048
WORKDIR=/tmp

k=23
m=13
psc_off=true
ps_skip=false
thr=0.7
ref_ind=/fs/cbcb-lab/rob/students/noor/Atacseq/piscem_analysis/hg38_ind_k${k}/hg38_ind_k${k}
piscem_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/piscem_analysis
out_start=$1
barcode_file=$2
read1_file=$3
read2_file=$4
base_filename=10k_pbmc_ATACv2_nextgem_Chromium_Controller_S3_L001
pesc_sc_atac=/fs/cbcb-lab/rob/students/noor/Atacseq/piscem_noor/piscem-cpp/build
out_dir=$piscem_dir/${out_start}_k${k}_psc_off=${psc_off}_ps_skip=${ps_skip}_thr=${thr}
mkdir -p $out_dir
#mkdir -p "$WORKDIR" && cd "$WORKDIR" || exit -1

if [[ -n $4 ]]
    then
        /usr/bin/time -o $out_dir/time_align.out $pesc_sc_atac/pesc-sc-atac --index $ref_ind \
        --read1 $read1_file \
        --read2 $read2_file \
        --barcode $barcode_file \
        --output $out_dir \
        --psc_off $psc_off \
        --ps_skip $ps_skip \
        --thr $thr
else
        /usr/bin/time -o $out_dir/time_align.out $pesc_sc_atac/pesc-sc-atac --index $ref_ind \
        --read1 $read1_file \
        --barcode $barcode_file \
        --output $out_dir \
        --psc_off $psc_off \
        --ps_skip $ps_skip \
        --thr $thr
fi

## Example
## sbatch helper_scripts/align_piscem.sh samp1285 sub_seq/barcode_subsamp.fq sub_seq/read1_subsamp.fq sub_seq/read2_subsamp.fq