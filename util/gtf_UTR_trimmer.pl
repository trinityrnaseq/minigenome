#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use GTF_utils;
use Fasta_reader;

my $usage = "\n\tusage: $0 file.gtf genome.fa > file.no_UTRs.gtf\n\n";

my $gtf_filename = $ARGV[0] or die $usage;
my $genome_fa = $ARGV[1] or die $usage;

main: {

    my $gene_obj_indexer = {};

    my $contig_to_genes = &GTF_utils::index_GTF_gene_objs($gtf_filename, $gene_obj_indexer);


    my $fasta_reader = new Fasta_reader($genome_fa);
    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    
    foreach my $contig (keys %$contig_to_genes) {

        my $contig_seq = $seqs{$contig};
        
        my @gene_ids = @{$contig_to_genes->{$contig}};

        foreach my $gene_id (@gene_ids) {

            my $gene_obj = $gene_obj_indexer->{$gene_id};
            
            $gene_obj->trim_UTRs();
            foreach my $isoform ($gene_obj->get_additional_isoforms()) {
                $isoform->trim_UTRs();
            }
            
            print $gene_obj->to_GTF_format(\$contig_seq);
            
        }
        
    }
    
    exit(0);

}
