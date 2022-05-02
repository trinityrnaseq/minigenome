#!/usr/bin/env python
# encoding: utf-8

#from __future__ import (absolute_import, division,
#                        print_function, unicode_literals)

import os, sys, re
import logging
import argparse
import pandas as pd

if sys.version_info[0] != 3:
    print("This script requires Python 3")
    exit(1)
            

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)



def main():

    parser = argparse.ArgumentParser(description="identify those transcript isoforms that are lacking genome representation", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--valid_stats", type=str, required=True, help="pasa_lite.align_stats.tsv")

    parser.add_argument("--genes_want", type=str, required=True, help="genes wanted")
    
    parser.add_argument("--isoform_expression_matrix",type=str, required=True, help="isoform expression matrix (TPM)")

    parser.add_argument("--min_sum_tpm", default=1.0, required=False, help="min sum TPM for a representative invalid isoform")
    

    args = parser.parse_args()

    
    valid_stats_filename = args.valid_stats
    genes_want_filename = args.genes_want
    isoform_expression_matrix_filename = args.isoform_expression_matrix
    min_sum_tpm = args.min_sum_tpm
    
        
    logger.info("parsing {}".format(isoform_expression_matrix_filename))
    expr_matrix = pd.read_csv(isoform_expression_matrix_filename, sep="\t", index_col=0)
    expr_sums = expr_matrix.sum(axis=1).reset_index()
    expr_sums.columns=['cdna_acc', 'sum_tpm']
    
    logger.info("parsing {}".format(valid_stats_filename))
    valid_stats = pd.read_csv(valid_stats_filename, sep="\t")


    logger.info("adding isoform expression to isoform stats")
    valid_stats = valid_stats.merge(expr_sums, how='left', on='cdna_acc')

    # clear out
    expr_matrix = None
    expr_sums = None
    

    logger.info("extracting gene ids from isoforms")
    def strip_isoform_val(iso_string):
        gene_string = re.sub("_i\\d+$", "", iso_string)
        return gene_string

    valid_stats['gene'] = valid_stats['cdna_acc'].apply(strip_isoform_val)
    logger.info(valid_stats.shape)

    logger.info("parsing {}".format(genes_want_filename))
    genes_want = pd.read_csv(genes_want_filename)
    logger.info(genes_want.shape)

    logger.info("filtering to retain only genes of interest")
    genes_want_stats = valid_stats.loc[ valid_stats['gene'].isin(genes_want['acc']) ]
    logger.info(genes_want_stats.shape)

    logger.info("filtering to retain those isoforms lacking any valid alignment")
    genes_want_stats = genes_want_stats.groupby('cdna_acc').filter(lambda x:  not x['valid_status'].isin(['VALID']).any() )
    logger.info(genes_want_stats.shape)


    
    # select a representative isoform (one with highest cumulative expression)
    logger.info("selecting representative isoform with highest cumulative expression")
    genes_want_stats = genes_want_stats.groupby('cdna_acc').apply(lambda x: x.sort_values('sum_tpm', ascending=False).head(1) )
    logger.info(genes_want_stats.shape)

    logger.info("pruning insufficiently expressed entries")
    genes_want_stats = genes_want_stats[gene_want_stats['sum_tpm'] >= min_sum_tpm ]
    logger.info(genes_want_stats.shape)
    
    logger.info("done, writing final output file of the missing isoforms")
    genes_want_stats.to_csv("invalid_isoforms_want.tsv", sep="\t", index=False)

    sys.exit(0)
    
 
####################
 
if __name__ == "__main__":
    main()
