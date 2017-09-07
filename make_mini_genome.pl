#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Nuc_translator;
use Overlap_piler;
use Data::Dumper;

my $max_intron_length = 5000;
my $min_intron_length = 200;
my $genome_flank_size = 1000;

my $out_prefix = "minigenome";


my $usage = <<__EOUSAGE__;

###############################################################################################
#
#  Required:
#
#  --gtf <string>                   genome annotation in gtf format 
#
#  --genome_fa <string>             genome sequence in fasta format
#
# Optional:
#
#  --shrink_introns
#
#  --min_intron_length <int>        default: $min_intron_length (only when --shrink_introns_used)
#  --max_intron_length <int>        default: $max_intron_length  (only when --shrink_introns used)
#
#  --genome_flank <int>             amt. of genomic sequence to extract flanking each gene (default: $genome_flank_size)
#
#  --out_prefix <string>            output prefix for output files (gtf and fasta) default: ${out_prefix}
#
#  --accs_restrict_file <string>    file containing lists of gene or transcript identifiers to restrict to
#
###############################################################################################


__EOUSAGE__

    ;

my $help_flag;

my $gtf_file;
my $genome_fasta_file;
my $shrink_introns_flag = 0;
my $accs_restrict_file = "";

&GetOptions ( 'h' => \$help_flag,
              
              'gtf=s' => \$gtf_file,
              'genome_fa=s' => \$genome_fasta_file,

              'shrink_introns' => \$shrink_introns_flag,
              'min_intron_length=i' => \$min_intron_length,
              'max_intron_length=i' => \$max_intron_length,
              'genome_flank=i' => \$genome_flank_size,
              
              'out_prefix=s' => \$out_prefix,

              'accs_restrict_file=s' => \$accs_restrict_file,
              
              
    );


if ($help_flag) {
    die $usage;
}

unless ($gtf_file && $genome_fasta_file) {
    die $usage;
}

main: {

    my %RESTRICT;

    if ($accs_restrict_file) {
        open(my $fh, $accs_restrict_file) or die $!;
        while (<$fh>) {
            chomp;
            my @x = split(/\s+/);
            foreach my $acc (@x) {
                $RESTRICT{$acc} = 1;
            }
        }
        close $fh;
    }
    
    
    my %gene_to_gtf = &extract_gene_gtfs($gtf_file);


    my @gene_info_structs = &make_gene_info_structs(%gene_to_gtf);
    
    open (my $out_genome_ofh, ">$out_prefix.fa.tmp") or die "Error, cannot write to $out_prefix.fa.tmp";
    open (my $out_gtf_ofh, ">$out_prefix.gtf.tmp") or die "Error, cannot write to $out_prefix.gtf.tmp";

    @gene_info_structs = sort {$a->{chr} cmp $b->{chr}
                               ||
                                   $a->{lend} <=> $b->{lend} } @gene_info_structs;
    

    my $supercontig = "";

    my $counter = 0;
    
    foreach my $gene_struct (@gene_info_structs) {

        $counter++; if ($counter > 1000) { last; }
                                
        my $gene_gtf = $gene_struct->{gtf};


        if (%RESTRICT) {

            unless (&includes_restricted_id($gene_gtf, \%RESTRICT)) {
                next;
            }
        }
                
        my ($gene_supercontig_gtf, $gene_sequence_region) = &get_gene_contig_gtf($gene_gtf, $genome_fasta_file);
        
        ($gene_supercontig_gtf, $gene_sequence_region) = &shrink_introns($gene_supercontig_gtf, $gene_sequence_region, $max_intron_length);
        
        if ($supercontig) {
            $supercontig .= ("N" x 1000);
        }
        
        $gene_supercontig_gtf = &adjust_gtf_coordinates($gene_supercontig_gtf, length($supercontig));

        $supercontig .= $gene_sequence_region;
        
        my $out_gtf = &set_gtf_scaffold_name("minigenome", $gene_supercontig_gtf);

        #$out_gtf = &include_gene_name_in_gene_id($out_gtf);
        
        print $out_gtf_ofh $out_gtf;
        
    }

    $supercontig =~ s/(\S{60})/$1\n/g;
    
    print $out_genome_ofh ">minigenome\n$supercontig\n";
    
    print STDERR "Done.\n";
    
    close $out_genome_ofh;
    close $out_gtf_ofh;

    if (! -s "$out_prefix.fa.tmp") {
        die "Error, no fusion contigs written";
    }
    else {
        rename("$out_prefix.fa.tmp", "$out_prefix.fa");
        rename("$out_prefix.gtf.tmp", "$out_prefix.gtf");
    }
    
    exit(0);
    

}

####
sub shrink_introns {
    my ($gene_gtf, $gene_seq_region, $max_intron_length) = @_;

    my @gtf_structs;
    my @gtf_lines = split(/\n/, $gene_gtf);
    foreach my $gtf_line (@gtf_lines) {
        my @x = split(/\t/, $gtf_line);
        push (@gtf_structs, \@x);
    }
    
    @gtf_structs = sort {$a->[3]<=>$b->[3]} @gtf_structs;

    
    ## get exon piles
    my $overlap_piler = new Overlap_piler();
    foreach my $gtf_row_aref (@gtf_structs) {
        
        my $exon_lend = $gtf_row_aref->[3];
        my $exon_rend = $gtf_row_aref->[4];
        $overlap_piler->add_coordSet($gtf_row_aref, $exon_lend, $exon_rend);
    }

    my @piles = $overlap_piler->build_clusters();
    
    my @pile_structs;
    foreach my $pile (@piles) {
        
        my @all_coords;
        foreach my $gtf_row_aref (@$pile) {
            my $lend = $gtf_row_aref->[3];
            my $rend = $gtf_row_aref->[4];
            push (@all_coords, $lend, $rend);
        }
        @all_coords = sort {$a<=>$b} @all_coords;

        my $pile_lend = shift @all_coords;
        my $pile_rend = pop @all_coords;
        
        my $pile_struct = { pile => $pile,
                            pile_lend => $pile_lend,
                            pile_rend => $pile_rend,
                            
                            pile_length => $pile_rend - $pile_lend + 1,

                    
                            new_pile_lend => $pile_lend,
                            new_pile_rend => $pile_rend,
        };
        push (@pile_structs, $pile_struct);
    }
   

    @pile_structs = sort { $a->{pile_lend} <=> $b->{pile_lend} } @pile_structs;

    ## set new pile bounds based on max intron length
    for (my $i = 1; $i <= $#pile_structs; $i++) {
        my $prev_pile_struct = $pile_structs[$i-1];
        my $curr_pile_struct = $pile_structs[$i];

        my $intron_length = $curr_pile_struct->{pile_lend} - $prev_pile_struct->{pile_rend} - 1;
        if ($intron_length > $max_intron_length) {
            $intron_length = $min_intron_length + int(rand($max_intron_length - $min_intron_length));
        }
        # set new relative coordinates for exon based on adjusted intron length
        $curr_pile_struct->{new_pile_lend} = $prev_pile_struct->{new_pile_rend} + $intron_length + 1;
        $curr_pile_struct->{new_pile_rend} = $curr_pile_struct->{new_pile_lend} + $curr_pile_struct->{pile_length} - 1;
        
    }
    
    ## adjust gtf exon coordinates
    
    my $gtf_adj = "";
    my $gene_seq_adj = "";
    
    my $prev_old_pile_rend = 0;
    my $prev_new_pile_rend = 0;
    
    foreach my $pile_struct (@pile_structs) {
        
        my $old_pile_lend = $pile_struct->{pile_lend};
        my $new_pile_lend = $pile_struct->{new_pile_lend};
        my $pile_length = $pile_struct->{pile_length};
        
        my $delta = $old_pile_lend - $new_pile_lend;
        
        ## add intron
        my $original_intron_len = $old_pile_lend - $prev_old_pile_rend -1;
        my $adj_intron_len = $new_pile_lend - $prev_new_pile_rend - 1;
        if ($adj_intron_len > $max_intron_length) {
            die "Error, intron length ($adj_intron_len) exceeds max of $max_intron_length";
        }
        my $intron_seq = "";
        if ($prev_old_pile_rend == 0 || $original_intron_len < $max_intron_length) {
            # use existing intron
            unless ($original_intron_len == $adj_intron_len) {
                die "Error, original short intron length is not equivalent to the adjusted intron length setting";
            }
            $intron_seq = substr($gene_seq_region, $prev_old_pile_rend, $original_intron_len);
        }
        else {
            ## split the difference
            my $left_intron_size = int($adj_intron_len/2);
            my $right_intron_size = $adj_intron_len - $left_intron_size;
            $left_intron_size -= 5;
            $right_intron_size -= 5; # add 10 Ns at center
            $intron_seq = substr($gene_seq_region, $prev_old_pile_rend, $left_intron_size);
            $intron_seq .= 'N' x 10;
            $intron_seq .= substr($gene_seq_region, $old_pile_lend - 1 - $right_intron_size, $right_intron_size);
            
            if (length($intron_seq) != $adj_intron_len) {
                die "Error, intron length is off: " . length($intron_seq) . " vs. $adj_intron_len (adjusted)";
            }
        }
        $gene_seq_adj .= $intron_seq;
        
        foreach my $gtf_row_aref (@{$pile_struct->{pile}}) {
            
            $gtf_row_aref->[3] -= $delta;
            $gtf_row_aref->[4] -= $delta;

            $gtf_adj .= join("\t", @$gtf_row_aref) . "\n";
        }

        my $pile_seq = substr($gene_seq_region, $old_pile_lend -1, $pile_length);
        $gene_seq_adj .= $pile_seq;

        $prev_old_pile_rend = $pile_struct->{pile_rend};
        $prev_new_pile_rend = $pile_struct->{new_pile_rend};
    }
    
    ## tack on end of sequence
    $gene_seq_adj .= substr($gene_seq_region, $prev_old_pile_rend);
    
    return($gtf_adj, $gene_seq_adj);
}


    
####
sub set_gtf_scaffold_name {
    my ($scaffold_name, $gtf_text) = @_;

    my $new_gtf = "";
    
    foreach my $line (split(/\n/, $gtf_text)) {
        
        my @x = split(/\t/, $line);
        $x[0] = $scaffold_name;
        
        #$x[8] =~ s/transcript_id \"/transcript_id \"$scaffold_name\^/;
        #$x[8] =~ s/gene_id \"/gene_id \"$scaffold_name\^/;
        
        $new_gtf .= join("\t", @x) . "\n";
    }
    
    return($new_gtf);
}

####
sub get_gene_contig_gtf {
    my ($gene_gtf, $genome_fasta_file) = @_;

    
    my ($gene_chr, $gene_lend, $gene_rend, $gene_orient, $revised_gene_gtf) = &get_gene_span_info($gene_gtf);
    
    $gene_gtf = $revised_gene_gtf;

    # print STDERR "INFO: gene_chr: $gene_chr, lend: $gene_lend, rend: $gene_rend\n";

    #print STDERR "\n\nGENE_GTF:\n$gene_gtf\n\n";
    
    
    my $seq_region = &get_genomic_region_sequence($genome_fasta_file,
                                                  $gene_chr, 
                                                  $gene_lend - $genome_flank_size, 
                                                  $gene_rend + $genome_flank_size,
                                                  $gene_orient);

    my $gene_contig_gtf = &transform_gtf_coordinates($gene_lend - $genome_flank_size,
                                                     $gene_gtf,
                                                     length($seq_region),
                                                     $gene_orient);
    

    # print STDERR "GENE_contig_gtf:\n$gene_contig_gtf\n\n";

    return($gene_contig_gtf, $seq_region);
    
    
}


#####
sub get_genomic_region_sequence {
    my ($fasta_file, $chr, $lend, $rend, $orient) = @_;

    my $cmd = "samtools faidx $fasta_file $chr:$lend-$rend";
    my $seq = `$cmd`;
    if ($?) {
        die "Error, cmd: $cmd died with ret $?";
    }
    my $header;
    ($header, $seq) = split(/\n/, $seq, 2);
    $seq =~ s/\s//g;
    
    my $seq_len = $rend - $lend + 1;
    if (length($seq) != $seq_len) {
        die "Error, didn't extract required sequence from $fasta_file, $chr, $lend, $rend, instead got seq of length " . length($seq);
    }
    
    return($seq);
}



####
sub extract_gene_gtfs {
    my ($gtf_file, $gene_want_href) = @_;

    my %gene_to_gtf;

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next;}
        unless (/\w/) { next; }
        my $line = $_;
        
        my @x = split(/\t/, $line);
        my $feat_type = $x[2];
        unless ($feat_type =~ /^(exon|CDS)$/) { next; } # only exon and CDS records

        my $gene_id = "";
        if (/gene_id \"([^\"]+)\"/) {
            $gene_id = $1;
        }
        
        my $chr = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        
        my $orig_info = "$chr,$gene_id,$lend,$rend,$orient";

        ## replace gene_id with gene_name
        if ($line =~ /gene_name \"([^\"]+)\"/) {
            my $gene_name = $1;
            $line =~ s/gene_id \"\Q$gene_id\E\"/gene_id \"$gene_name\"/;
            $gene_id = $gene_name;

            $line =~ s/transcript_id \"/transcript_id \"$gene_name^/;
            
        }
        
        $line .= " orig_coord_info \"$orig_info\";\n";
        
        $gene_to_gtf{$gene_id} .= $line;
                
    }
    close $fh;

    return(%gene_to_gtf);
}


####
sub get_gene_span_info {
    my ($gene_gtf_text) = @_;

    my ($chr, $min_lend, $max_rend, $orient);

    my $revised_gene_gtf_text = "";
    my $template_line = "";

    my @gtf_lines = split(/\n/, $gene_gtf_text);
    foreach my $line (@gtf_lines) {
        my @x = split(/\t/, $line);
        my $scaffold = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        my $strand = $x[6];
        if (defined $chr) {
            ## check to ensure the rest of the info matches up
            ## if discrepancy, use first found.
            if ($chr ne $scaffold) {
                print STDERR "Error, chr discrepancy in gtf info [chr: [$chr] vs. [$scaffold] ]for $line\n\ttemplate: $template_line\n\n";
                next; 
            }
            if ($orient ne $strand) {
                print STDERR "Error, strand conflict [$orient] vs. [$strand] in gtf info for $line\n\ttemplate: $template_line\n\n";
                next;
            }
            if ($lend < $min_lend) {
                $min_lend = $lend;
            }
            if ($rend > $max_rend) {
                $max_rend = $rend;
            }
            
        }
        else {
            ## init
            ($chr, $min_lend, $max_rend, $orient) = ($scaffold, $lend, $rend, $strand);
            $template_line = "$line";
        }
        $revised_gene_gtf_text .= "$line\n";
    }

    if ($max_rend - $min_lend > 100e6) {
        confess "Error - no gene spans 100M bases in length.... likely problem";
    }
    
    return($chr, $min_lend, $max_rend, $orient, $revised_gene_gtf_text);
}


####
sub transform_gtf_coordinates {
    my ($left_reference_pos, $gene_gtf, $seq_length, $gene_orient) = @_;

    my $new_gtf = "";
    
    foreach my $line (split(/\n/, $gene_gtf)) {
        
        my @fields = split(/\t/, $line);
        my ($lend, $rend) = ($fields[3], $fields[4]);
        
        $lend = $lend - $left_reference_pos + 1;
        $rend = $rend - $left_reference_pos + 1;
        
        #if ($gene_orient eq '-') {
            # revcomp the coordinates
        #    $lend = $seq_length - $lend + 1;
        #    $rend = $seq_length - $rend + 1;
        #    ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        #}

        $fields[3] = $lend;
        $fields[4] = $rend;
        #$fields[6] = '+';
        
        $new_gtf .= join("\t", @fields) . "\n";
    }
    
    return($new_gtf);

}

####
sub adjust_gtf_coordinates {
    my ($gtf, $adjustment) = @_;

    my $new_gtf = "";
   
    foreach my $line (split(/\n/, $gtf)) {

        my @fields = split(/\t/, $line);
        $fields[3] += $adjustment;
        $fields[4] += $adjustment;

        $new_gtf .= join("\t", @fields) . "\n";
    }

    return($new_gtf);
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    return;
}

####
sub make_gene_info_structs {
    my (%gene_to_gtf) = @_;

    my @gene_structs;
    
    foreach my $gene (keys %gene_to_gtf) {
        my $gtf = $gene_to_gtf{$gene};

        my ($chr, $min_lend, $max_rend, $orient, $revised_gene_gtf_text) = &get_gene_span_info($gtf);

        my $struct = { chr => $chr,
                       lend => $min_lend,
                       rend => $max_rend,
                       orient => $orient,
                       gtf => $revised_gene_gtf_text,
                       gene => $gene,
        };

        push (@gene_structs, $struct);
    }

    return(@gene_structs);
}


####
sub includes_restricted_id {
    my ($gene_gtf, $restrict_href) = @_;

    $gene_gtf =~ /gene_id \"([^\"]+)\"/ or die "Error, cannot extract gene_id from $gene_gtf";
    my $gene_id = $1;

    if ($restrict_href->{$gene_id}) {
        return(1);
    }

    if ($gene_gtf =~ /gene_name \"([^\"]+)\"/) {
        my $gene_name = $1;
        if ($restrict_href->{$gene_name}) {
            return(1);
        }
    }
    
    while ($gene_gtf =~ /transcript_id \"([^\"]+)\"/g) {
        my $tx_id = $1;
        if ($restrict_href->{$tx_id}) {
            return(1);
        }
    }

    return(0);
}


####
sub include_gene_name_in_gene_id {
    my ($gene_gtf) = @_;

    
    if ($gene_gtf =~ /gene_name \"([^\"]+)\"/) {
        my $gene_name = $1;

        $gene_gtf =~ /gene_id \"([^\"]+)\"/ or die "Error, not extracting gene_id from $gene_gtf";
        my $gene_id = $1;

        #print "Gene_name: $gene_name, gene_id: $gene_id\n";

        #print "Before: $gene_gtf\n";
        #$gene_gtf =~ s/gene_id \"${gene_id}/gene_id \"${gene_name}\^${gene_id}/g or die "Error, no name update for $gene_gtf";
        $gene_gtf =~ s/gene_id \"\Q${gene_id}\E\"/gene_id \"$gene_name^$gene_id\"/g or die "Error, no name update for $gene_gtf";
        #print "After: $gene_gtf\n\n";
        
    }
    
            

    return($gene_gtf);

}

