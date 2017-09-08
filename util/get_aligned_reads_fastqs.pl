#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Fastq_reader;
use Process_cmd;
use Data::Dumper;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);  


my $usage = <<__EOUSAGE__;

######################################################################################
#
#  --bam <string>               aligned reads bam file
#
#  --left_fq <string>           /1 fastq file (or SE reads)
#
#  --right_fq <string>          /2 fastq file (optional, required for PE reads)
#
#  --output_prefix <string>     output prefix
#
######################################################################################

__EOUSAGE__

    ;


my $help_flag;

my $bam_file;
my $left_fq;
my $right_fq;
my $output_prefix;


&GetOptions( 'help|h' => \$help_flag,

             'bam=s' => \$bam_file,
             
             'left_fq=s' => \$left_fq,
             'right_fq=s' => \$right_fq,
             'output_prefix=s' => \$output_prefix,
    );

if ($help_flag) {
    die $usage;
}

unless ($bam_file && $left_fq && $output_prefix) {
    die $usage;
}


main: {

    my %aligned_read_names = &get_aligned_read_names($bam_file);
        
    &write_fastq_files($left_fq, "$output_prefix.aligned_reads_1.fq", \%aligned_read_names);
    
    if ($right_fq) {
        &write_fastq_files($right_fq, "$output_prefix.aligned_reads_2.fq", \%aligned_read_names);
    }
    
    print STDERR "\nDone.\n\n";
    
    exit(0);
    
}

####
sub get_aligned_read_names {
    my ($bam_file) = @_;

    my %aligned_reads;
    
    open(my $fh, "samtools view $bam_file | ") or die "Error, cannot read bam file: $bam_file";
    while (<$fh>) {
        my @x = split(/\t/);
        my $acc = $x[0];
        $acc =~ s/\/[12]$//;
        my $cigar = $x[5];
        if ($cigar =~ /M/) {
            $aligned_reads{$acc} = 1;
        }
    }
    close $fh;

    return(%aligned_reads);
}


####
sub write_fastq_files {
    my ($input_fastq_file, $output_fastq_file, $aligned_reads_href) = @_;
    
    my $fastq_reader = new Fastq_reader($input_fastq_file);

    open (my $ofh, ">$output_fastq_file") or die "Error, cannot write to $output_fastq_file";
    
    while (my $fq_record = $fastq_reader->next()) {
        
        my $core_read_name = $fq_record->get_core_read_name();
        #print "[$core_read_name]\n";
        
        if (my $fusion_name = $aligned_reads_href->{$core_read_name}) {
            
            my $record_text = $fq_record->get_fastq_record();
            chomp $record_text;
            my @lines = split(/\n/, $record_text);
            my ($_1, $_2, $_3, $_4) = @lines;
            print $ofh join("\n", ($_1, $_2, $_3, $_4)) . "\n";
        }
    }
    
    print STDERR "\nDone writing to $output_fastq_file\n\n";
    
    return;
}


####
sub append_reads_to_fusion {
    my ($fusion_name, $core_frag_name_to_fusion_name_href, $reads_href) = @_;

    foreach my $frag_name (keys %$reads_href) {
        $core_frag_name_to_fusion_name_href->{$frag_name} = $fusion_name ;
    }
    
    return;
}

####
sub parse_core_frag_names {
    my ($comma_delim_read_list_txt) = @_;

    my %core_frag_names;

    my @read_names = split(/,/, $comma_delim_read_list_txt);
    foreach my $read_name (@read_names) {
        $read_name =~ s|/[12]$||;
        $core_frag_names{$read_name} = 1;
    }

    
    return(%core_frag_names);
}


