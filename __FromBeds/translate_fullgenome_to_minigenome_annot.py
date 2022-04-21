#!/usr/bin/env python3

import sys, os, re
import subprocess
import logging
import argparse
import pandas as pd
import pyranges as pr
import csv

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    #parser = argparse.ArgumentParser(description="translate full genome bed to the minigenome bed",
     #                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    #parser.add_argument("--fullgenome_annot", type=str, required=True, help="input fullgenome bed or gtf file to translate coords")
    #parser.add_argument("--translation_intervals_bed", type=str, required=True, help="intervals bed file")

    #args = parser.parse_args()

    
    df = pd.read_csv("AmexT_v47-AmexG_v6.0-DD.no_seqs.chr10p.gtf", sep="\t", names=["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info"])
    pd_translation = pd.read_csv("minigenome.coord_translation.tsv", sep="\t")

    pr_df = pr.PyRanges(df)
    pr_translation = pr.PyRanges(pd_translation)
    
    pr_join = pr_df.join(pr_translation)

    pr_join.df.to_csv("AmexT_v47-AmexG_v6.0-DD.no_seqs.gtf.w_translation_table.tsv", sep="\t", index=False, quoting=csv.QUOTE_NONE)

    data = pd.read_csv("AmexT_v47-AmexG_v6.0-DD.no_seqs.gtf.w_translation_table.tsv", sep="\t")
    data['Start2'] = data['Start'] - data['Start_b'] + data['NewStart']
    data['End2'] = data['End'] - data['Start_b'] + data['NewStart']

    data['OrigStart'] = data['Start']
    data['OrigEnd'] = data['End']
    
    data['Start'] = data['Start2']
    data['End'] = data['End2']
    data = data[ ["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info", "NewStart", "OrigStart", "OrigEnd"] ]

    #data = data[ ["Chromosome", "Source", "Type", "Start", "End", "Something", "Start_b", "End_b", "NewStart", "NewEnd", "OrigStart", "OrigEnd", "Strand", "Dot", "Info"] ]
                 
    data.to_csv("AmexT_v47-AmexG_v6.0-DD.no_seqs.minigenome.gtf", sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)

    


    """
    minigenome = "minigenome.fa"
    fullgenome = "chr10p.fa"


    full_chrom_seq = extract_genome_seq(fullgenome, "chr10p")
    mini_chrom_seq = extract_genome_seq(minigenome, "chr10p") 

    
    for i, row in data.iterrows():
        # validate feature seq
        # Start,End seq should matchj OrigStart,OrigEnd

        orig_seq = full_chrom_seq[row['OrigStart']-1:row['OrigEnd']]
        mini_seq = mini_chrom_seq[row['Start']-1:row['End']]
        
        if orig_seq != mini_seq:
            print(f"row {row} !=  orig_seq: {orig_seq} vs. mini_seq: {mini_seq}")
        else:
            print("OK")
        


    """

    """

    pd_annot = pd.read_csv("AmexT_v47-AmexG_v6.0-DD.gtf.exons_n_CDS.bed", sep="\t", names=["Chromosome", "Start", "End"])
    pd_annot = pd_annot.drop_duplicates()
    pd_translation = pd.read_csv("minigenome.bed", sep="\t")

    chromosomes = pd_translation.Chromosome.unique()


    pd_final = None

    for chromosome in chromosomes:

        logger.info("-processing {}".format(chromosome))

        pd_annot_chr = pd_annot[pd_annot.Chromosome == chromosome]
        pd_translation_chr = pd_translation[pd_translation.Chromosome == chromosome]
        
        pr_annot_chr = pr.PyRanges(pd_annot_chr)
        pr_translation_chr = pr.PyRanges(pd_translation_chr)

        pr_annot_chr = pr_annot_chr.join(pr_translation_chr)

        pd_final = pd.concat([pd_final, pr_annot_chr.df])


    pd_final.to_csv("pd_final.csv", sep="\t")

    """
    
    #pr_annot = pr.PyRanges(pd_annot)
    #pr_annot = pr_annot.join(pr_translation)
    
    
    


    #other_bed_df = pd.read_csv("AmexG_MINIGENOME/pd_final.csv", sep="\t")

    #j = df.merge(other_bed_df, left_on=["Chromosome", "Lend", "Rend"], right_on=["Chromosome", "Start", "End"])

    

    #fullgenome_annot_file = args.fullgenome_annot
    #translation_intervals_bed_file = args.translation_intervals_bed



    """

    logger.info("parsing {}".format(fullgenome_annot_file))
    if re.search("\\.gtf$", fullgenome_annot_file):
        pr_annot = pr.read_gtf(fullgenome_annot_file)
    elif re.search("\\.bed$", fullgenome_annot_file):
        pr_annot = pr.read_bed(fullgenome_annot_file)
    else:
        raise RuntimeError("not recognizing {} as bed or gtf file".format(fullgenome_annot_file))

    
    translation_intervals = pd.read_csv(translation_intervals_bed, sep="\t")
    pr_translation_intervals = pr.PyRanges(translation_intervals)

    """



def extract_genome_seq(genome_fasta_file, chromosome):
    
    cmd = f"samtools faidx {genome_fasta_file} {chromosome}"
    logger.info(cmd)
    fa_record = subprocess.check_output(cmd, shell=True).decode()
    fa_record = fa_record.split("\n")
    fa_record = fa_record[1:]
    fa_record = "".join(fa_record)
    
    return fa_record



if __name__=='__main__':
    main()