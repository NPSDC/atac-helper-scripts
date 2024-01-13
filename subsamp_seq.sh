#!/bin/bash

fastq_dir=/fs/cbcb-lab/rob/students/noor/Atacseq/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fastqs
base_filename=10k_pbmc_ATACv2_nextgem_Chromium_Controller_S3_L001
seqkit_path=/fs/cbcb-lab/rob/students/noor/Atacseq/seqkit 

$seqkit_path sample -p 0.000005 $fastq_dir/${base_filename}_R1_001.fastq.gz -j 16  -o sub_seq/read1_subsamp.fq
$seqkit_path sample -p 0.000005 $fastq_dir/${base_filename}_R3_001.fastq.gz -j 16  -o sub_seq/read2_subsamp.fq
$seqkit_path sample -p 0.000005 $fastq_dir/${base_filename}_R2_001.fastq.gz -j 16  -o sub_seq/barcode_subsamp.fq

### We can have barcodes repeated ideally (each barcode corresponds to cell)
### Currently we sre not removing duplicated barcode
