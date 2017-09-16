#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Process_cmd;
    
my $usage = "\n\n\tusage: $0 minigenome.fa minigenome.gtf samples_file.txt\n\n";

my $minigenome_fa = $ARGV[0] or die $usage;
my $annotation_gtf = $ARGV[1] or die $usage;
my $samples_file = $ARGV[2] or die $usage;

main: {

    my @read_sets = &parse_samples_file($samples_file);
    
    ## align reads to the mini-genome using hisat2

    unless (-e "$minigenome_fa.hisat2.build.ok") {
    
        &process_cmd("hisat2_extract_splice_sites.py $annotation_gtf  > $annotation_gtf.ss");
        
        &process_cmd("hisat2_extract_exons.py $annotation_gtf > $annotation_gtf.exons");
        
        &process_cmd("hisat2-build --exon $annotation_gtf.exons --ss $annotation_gtf.ss  $minigenome_fa $minigenome_fa");

        &process_cmd("touch $minigenome_fa.hisat2.build.ok");
    }
    
    foreach my $read_set_aref (@read_sets) {
        my ($sample_id, $left_fq, $right_fq) = @$read_set_aref;

        
        &process_cmd("set -eof pipefail; hisat2 --dta -x $minigenome_fa -1 $left_fq -2 $right_fq | samtools view -Sb -F 4 | samtools sort -o $sample_id.hisat2.bam");

        ## extract the reads that map to the mini-genome:

        &process_cmd("$FindBin::Bin/get_aligned_reads_fastqs.pl --bam $sample_id.hisat2.bam --left_fq $left_fq --right_fq $right_fq --output_prefix $sample_id");
    }

    exit(0);
    
}



####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;
    
    open(my $fh, $samples_file) or die "Error, cannot open file $samples_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($cond, $rep, $fq_a, $fq_b) = @x;

        my $readset_name = "${cond}-${rep}";
        
        push (@samples, [$readset_name, $fq_a, $fq_b]);
    }
    close $fh;

    return (@samples);
}

