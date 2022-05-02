#!/usr/bin/env python3

import sys, os, re
import subprocess
import logging
import argparse
import pandas as pd
import pyranges as pr

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="build supplemental transcript contigs fasta",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--transcripts_fasta", type=str, required=True, help="fasta file containing all transcripts")
    parser.add_argument("--transcript_acc_list", type=str, required=True, help="file containing a list of tanscript identifiers to use (column name should be 'cdna_acc').")
    parser.add_argument("--Nspacer_length", type=int, default=10, help="spacer length between transcripts")
    parser.add_argument("--contig_name", type=str, default="supp_transcripts", help="contig name and prefix for output files: .genome.fa and .genome.bed")
    
    args = parser.parse_args()

    transcripts_fasta_filename = args.transcripts_fasta
    transcript_acc_list_file = args.transcript_acc_list
    Nspacer_length = args.Nspacer_length
    contig_name = args.contig_name


    transcript_info = pd.read_csv(transcript_acc_list_file, sep="\t")
    transcript_accs = transcript_info['cdna_acc'].tolist()
    transcript_accs = set(transcript_accs)

    transcript_seqs = get_transcript_seqs(transcript_accs, transcripts_fasta_filename)


    genome_fa_ofh = open(contig_name + ".genome.fa", "wt")
    genome_bed_ofh = open(contig_name + ".genome.bed", "wt")
    
    Nspacer = "N" * Nspacer_length

    contig_new_lend = 0
    contig_new_rend = 0

    
    print(f">{contig_name}", file=genome_fa_ofh)

    transcript_accs = sorted(list(transcript_accs))
    
    for transcript_acc in transcript_accs:
   
        transcript_seq = transcript_seqs[transcript_acc]
        transcript_len = len(transcript_seq)

        contig_new_lend += 1
        contig_new_rend = contig_new_lend + transcript_len - 1
            
        print("\t".join([contig_name, "supercontig", "cDNA_match", str(contig_new_lend), str(contig_new_rend),
                         ".", "+", ".",
                         f"ID={transcript_acc};Parent={transcript_acc};Target={transcript_acc} 1 {transcript_len}"]),
              file=genome_bed_ofh)

        # add short N spacer
        print(transcript_seq + Nspacer, file=genome_fa_ofh, end='')
            
        contig_new_lend = contig_new_rend + Nspacer_length
            

    print("", file=genome_fa_ofh) # newline seq separator
    
    
    logger.info("done")
        
    sys.exit(0)




def get_transcript_seqs(acc_set, fasta_filename):

    acc_set = acc_set.copy()
    
    transcript_seqs = dict()

    transcript_seq = ""
    acc = ""
    with open(fasta_filename, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            m = re.match(">(\S+)", line)
            if m:
                next_acc = m.group(1)

                if acc in acc_set:
                    transcript_seqs[acc] = transcript_seq
                    acc_set.remove(acc)
                    
                transcript_seq = ""
                acc = next_acc
            else:
                transcript_seq += line

    # get last one
    if acc in acc_set:
        transcript_seqs[acc] = transcript_seq
        acc_set.remove(acc)

    if len(acc_set) == 0:
        logger.info("all transcript seqs captured")
    else:
        raise RuntimeError("Error, missing " + str(len(acc_set)) + " sequences from transcripts fasta file")
    
    return transcript_seqs



if __name__=='__main__':
    main()
