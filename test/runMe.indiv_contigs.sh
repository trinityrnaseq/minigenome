#!/bin/bash

set -evo pipefail

## make mini-genome
../make_mini_genome.pl  --gtf data/annotation.gtf --genome_fa data/genome.fa --shrink_introns --accs_restrict_file genes.want.txt --single_contig_per_gene --out_prefix mini_indiv

## align reads to the mini-genome using hisat2

hisat2_extract_splice_sites.py mini_indiv.gtf  > mini_indiv.gtf.ss

hisat2_extract_exons.py  mini_indiv.gtf > mini_indiv.gtf.exons

hisat2-build --exon ./mini_indiv.gtf.exons --ss ./mini_indiv.gtf.ss  ./mini_indiv.fa ./mini_indiv.fa

hisat2 --dta -x mini_indiv.fa -1 data/reads_1.fq -2 data/reads_2.fq | samtools view -Sb -F 4 | samtools sort -o alignments.mini_indiv.hisat2.bam

## extract the reads that map to the mini-genome:

../util/get_aligned_reads_fastqs.pl  --bam alignments.mini_indiv.hisat2.bam --left_fq data/reads_1.fq --right_fq data/reads_2.fq --output_prefix mini_indiv

# done. :-)



