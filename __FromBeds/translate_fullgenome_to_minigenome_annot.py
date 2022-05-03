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

    parser = argparse.ArgumentParser(description="translate full genome gtf|gff3 to the minigenome gtf|gff3",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument("--fullgenome_annot", type=str, required=True, help="input fullgenome gtf or gff3 file to translate coords")
    parser.add_argument("--translation_intervals", type=str, required=True, help="translation intervals tsv file")
    parser.add_argument("--output_gtf", type=str, required=True, help="output gtf filename")


    args = parser.parse_args()

    fullgenome_annot_filename = args.fullgenome_annot
    translation_intervals_filename = args.translation_intervals
    output_gtf = args.output_gtf
    

    logger.info("parsing {}".format(fullgenome_annot_filename))

    df = pd.read_csv(fullgenome_annot_filename, sep="\t", names=["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info"])


    logger.info("parsing {}".format(translation_intervals_filename))
    pd_translation = pd.read_csv(translation_intervals_filename, sep="\t")


    logger.info("joining based on feature coordinate overlaps")
    pr_df = pr.PyRanges(df)
    pr_translation = pr.PyRanges(pd_translation)    
    pr_join = pr_df.join(pr_translation)

    logger.info("assigning adjusted coordinates")
    data = pr_join.df.copy()
    data['Start2'] = data['Start'] - data['Start_b'] + data['NewStart']
    data['End2'] = data['End'] - data['Start_b'] + data['NewStart']

    data['OrigStart'] = data['Start']
    data['OrigEnd'] = data['End']

    def add_orig_coords(row):
        return row['Info'] + f" OrigCoords \"{row['OrigStart']}-{row['OrigEnd']}\" "
    
    data['Info'] = data.apply(add_orig_coords, axis=1)
    
        
    
    data['Start'] = data['Start2']
    data['End'] = data['End2']
    data = data[ ["Chromosome", "Source", "Type", "Start", "End", "Something", "Strand", "Dot", "Info"] ]

    logger.info("writing output gtf: {}".format(output_gtf))
    data.to_csv("minigenome.gtf", sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='\\')


    sys.exit(0)


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
