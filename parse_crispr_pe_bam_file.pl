#! /usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $help;
my $label;
my $input_bam_file;
my $output_pdf_file;
my $debug;
my $ref_seq_file;
my $min_overlap = 10;
my $min_length = 80;
my $desc = qq{parse_crispr_pe_bam_file.pl [options] --input in.bam --output out.pdf

Description:
This script reads a bam file with PE reads and check that the reads overlap and that they match on
the overlap. All pairs that pass these check are used to look for insertions and deletions.

The output is a set of PDF, PNG, SVG and TXT files with the result of the analysis, mostly different
plots showing the different deletions and the list of the most common deletions.

Options:
--help: shows this help
--debug: output debugging information
--ref_seq: fasta reference sequence file

};

GetOptions(
    "help"  => \$help,
    "debug"  => \$debug,
    "label=s" => \$label,
    "input_file=s" => \$input_bam_file,
    "output_file=s" => \$output_pdf_file,
    "ref_seq=s" => \$ref_seq_file,
    );

if ($help or !$input_bam_file or !$output_pdf_file) {
    print $desc;
    exit();
}

if (!$label) {
    $label = $input_bam_file;
}

my ($data_file, $stats) = parse_bam_file($input_bam_file);

my $total = 0;
foreach my $key (sort keys %$stats) {
    my $value = $stats->{$key};
    print "$key: $value\n";
    $total += $value;
}
print "TOTAL: $total\n";

my $wt = ($stats->{"OK:WILD-TYPE"} or 0);

my $top_sequences = [];
if ($ref_seq_file) {
    $top_sequences = get_top_sequences($ref_seq_file, $data_file);
}

my $R_script = write_R_script($label, $data_file, $top_sequences, $wt, $output_pdf_file);

print qx"Rscript $R_script";

exit(0);


=head2 get_top_sequences

  Arg[1]        : string $ref_seq_filename
  Arg[2]        : string $data_filename
  Example       : my $top_sequences = get_top_sequences($ref_seq_file, $$data_file);
  Description   : Reads from the $data_file the most common deletions and insertions (up to 10) and
                  aligns them to the wild-type sequence.
  Returns       : arrayref of strings
  Exceptions    : Dies if any filename is not found.

=cut

sub get_top_sequences {
    my ($ref_seq_file, $data_file) = @_;
    my $top_sequences = [];

    ## ------------------------------------------------------------------------------
    ## Reads the fasta sequence
    ## ------------------------------------------------------------------------------
    open(REF, $ref_seq_file) or die;
    my $ref_seq;
    while (<REF>) {
        chomp;
        next if (/^>/);
        $ref_seq .= $_;
    }


    ## ------------------------------------------------------------------------------
    ## Reads and extract stats on most common deletions and insertions:
    ## ------------------------------------------------------------------------------
    my @del_lines = qx"more $data_file | awk '\$2 == \"DEL\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn | head -n 10";
    my @ins_lines = qx"more $data_file | awk '\$2 == \"INS\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn | head -n 10";


    ## ------------------------------------------------------------------------------
    ## Calculates the region of the REF sequence to display
    ## ------------------------------------------------------------------------------
    my $min_from;
    my $max_to;
    my $longest_seq = 0;
    foreach my $this_line (@del_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        $min_from = $from if (!$min_from or $from < $min_from);
        $max_to = $from+$del_length if (!$max_to or $from+$del_length > $max_to);
        $longest_seq = length($seq) if ($longest_seq < length($seq));
    }
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        $min_from = $from-1 if (!$min_from or $from-1 < $min_from);
        $max_to = $from if (!$max_to or $from > $max_to);
        $longest_seq = length($seq)+2 if ($longest_seq < length($seq)+2);
    }
    $min_from -= 21;
    $min_from = 0 if ($min_from < 0);
    $max_to += 19;

    ## ------------------------------------------------------------------------------
    ## Sets the output format (using whitespaces to be nicely printed in R afterewards)
    ## ------------------------------------------------------------------------------
    my $format = "\%-".($max_to-$min_from)."s \%7s \%4s \%3s \%6s \%-${longest_seq}s";


    ## ------------------------------------------------------------------------------
    ## Header and WT sequence
    ## ------------------------------------------------------------------------------
    my $header = sprintf($format, "Sequence" , "Num", "TYPE", "L", "POS", "Diff");
    my $wt_sequence = sprintf($format, substr($ref_seq, $min_from, $max_to-$min_from), $wt, "W-T", 0, "NA", "");
    $top_sequences = [$header, "", $wt_sequence, ""];


    ## ------------------------------------------------------------------------------
    ## Most common deletions
    ## ------------------------------------------------------------------------------
    foreach my $this_line (@del_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $from-$min_from-1) . '-'x$del_length . substr($ref_seq, ($from+$del_length-1), $max_to - ($from+$del_length-1));
        push(@$top_sequences, sprintf($format, $aligned_seq, $num, "DEL", $del_length, $from, $seq));
    }


    ## ------------------------------------------------------------------------------
    ## Separation line
    ## ------------------------------------------------------------------------------
    push(@$top_sequences, "");


    ## ------------------------------------------------------------------------------
    ## Most common insertions
    ## ------------------------------------------------------------------------------
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $max_to-$min_from);
        substr($aligned_seq, $from-$min_from-2, 2, "><");
        push(@$top_sequences, sprintf($format, $aligned_seq, $num, "INS", $ins_length, $from, ">$seq<"));
    }

    # Print this on the standard output
    print "\n", join("\n", @$top_sequences), "\n";

    # Return the lines (as an arrayref)
    return $top_sequences;
}


=head2 parse_bam_file

  Arg[1]        : string $bam_filename
  Example       : my ($data_filename, $stats) = parse_bam_file($bam_file);
  Description   : Extracts from the $bam_file all the PE reads that match all the criteria (aligned,
                  > min_length, good overlap, no indel or mismatch before or after the overlap,
                  perfect match between both reads and just a single insertion or deletion (but no
                  mismatches) w.r.t. the ref sequence in the overlap
  Returns       : array of $data_filename and $stats. These are:
                  - string with the filename where the resulting insertion and deletions have been
                    stored
                  - hashref of keys (event) and values (number of these events)
  Exceptions    : Dies if it cannot open the file with samtools

=cut

sub parse_bam_file {
    my ($bam_file) = @_;
    my $stats;

    open(SAM, "samtools view $bam_file |") or die;

    my $data_file = "$bam_file.data.txt";
    open(DATA, ">$data_file");
    print DATA join("\t", "read", "event", "event_length", "from", "to", "midpoint", "seq"), "\n";

    while (<SAM>) {

        ## ------------------------------------------------------------------------------
        ## Read both lines (assuming unsorted BAM file)
        ## ------------------------------------------------------------------------------
        chomp;
        my ($qname1, $flag1, $rname1, $pos1, $mapq1, $cigar1, $rnext1, $pnext1, $tlen1, $seq1, $qual1, @others1) = split("\t", $_);
        $_ = <SAM>;
        chomp;
        my ($qname2, $flag2, $rname2, $pos2, $mapq2, $cigar2, $rnext2, $pnext2, $tlen2, $seq2, $qual2, @others2) = split("\t", $_);
        my $md1 = (grep {/^MD:Z:/} @others1)[0];
        my $md2 = (grep {/^MD:Z:/} @others2)[0];
        my $end1 = $pos1 + length($seq1);
        my $end2 = $pos2 + length($seq2);
        if ($debug) {
          print join("\t", $qname1, $flag1, $rname1, $pos1, $mapq1, $cigar1, $rnext1, $pnext1, $tlen1), "\n";
          print join("\t", $qname2, $flag2, $rname2, $pos2, $mapq2, $cigar2, $rnext2, $pnext2, $tlen2), "\n";
        }


        ## ------------------------------------------------------------------------------
        ## Check that both lines refer to the same read
        ## ------------------------------------------------------------------------------
        if ($qname1 ne $qname2) {
            die "SAM file sorted or not for PE reads\n";
        }
        die "SAM file seems to have been sorted in some way. PE reads are not consecutive\n" if ($rname1 ne $rname2);


        ## ------------------------------------------------------------------------------
        ## Check that both reads align (i.e., they have a cigar string)
        ## ------------------------------------------------------------------------------
        if ($cigar1 eq "*" or $cigar2 eq "*") {
            $stats->{"01.no_cigar"}++;
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Check that both reads are long enough
        ## ------------------------------------------------------------------------------
        if (length($seq1) < $min_length or length($seq2) < $min_length) {
            $stats->{"02.short_read"}++;
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Check that both reads overlap and that the overlap is long enough
        ## ------------------------------------------------------------------------------
        my $overlap_start = $pos1>$pos2?$pos1:$pos2;
        my $overlap_end = $end1<$end2?$end1:$end2;
        if ($overlap_end < $overlap_start + $min_overlap) {
            $stats->{"03.no_overlap"}++;
            next;
        }
#        print "Overlap: $overlap_start-$overlap_end\n";


        ## ------------------------------------------------------------------------------
        ## Keeps a safe copy of both original sequences
        ## ------------------------------------------------------------------------------
        my $original_seq1 = $seq1;
        my $original_seq2 = $seq2;


        ## ------------------------------------------------------------------------------
        ## Clip the sequence 3' of the overlap (and check that there is no mismatch nor indel)
        ## ------------------------------------------------------------------------------
        if ($pos1 > $pos2) {
            # Clip R2
            my $diff_in_ref_bp = $pos1 - $pos2;
            my ($initial_match) = $cigar2 =~ /^(\d*)M/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if ($initial_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"05.indel_in_3prime"}++;
                next;
            }
            ($initial_match) = $md2 =~ /^MD:Z:(\d*)/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if ($initial_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"06.mismatch_in_3prime"}++;
                next;
            }
            substr($seq2, 0, $diff_in_ref_bp, "");
            print "1A $seq1\n1A $seq2\n" if ($debug);
        } elsif ($pos1 < $pos2) {
            # Clip R1
            my $diff_in_ref_bp = $pos2 - $pos1;
            my ($initial_match) = $cigar1 =~ /^(\d*)M/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if ($initial_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"05.indel_in_3prime"}++;
                next;
            }
            ($initial_match) = $md1 =~ /^MD:Z:(\d*)/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if ($initial_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"06.mismatch_in_3prime"}++;
                next;
            }
            substr($seq1, 0, $diff_in_ref_bp, "");
            print "1B $seq1\n1B $seq2\n" if ($debug);
        }


        ## ------------------------------------------------------------------------------
        ## Clip the sequence 5' of the overlap (and check that there is no mismatch nor indel)
        ## ------------------------------------------------------------------------------
        if ($end1 > $end2) {
            # Clip R1
            my $diff_in_ref_bp = $end1 - $end2;
            my ($last_match) = $cigar1 =~ /(\d*)M$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if ($last_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"07.indel_in_5prime"}++;
                next;
            }
            ($last_match) = $md1 =~ /(\d*)$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if ($last_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"08.mismatch_in_5prime"}++;
                next;
            }
            substr($seq1, -$diff_in_ref_bp, $diff_in_ref_bp, "");
            print "2A $seq1\n2A $seq2\n" if ($debug);
        } elsif ($end1 < $end2) {
            # Clip R2
            my $diff_in_ref_bp = $end2 - $end1;
            my ($last_match) = $cigar2 =~ /(\d*)M$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if ($last_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"07.indel_in_5prime"}++;
                next;
            }
            ($last_match) = $md2 =~ /(\d*)$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if ($last_match < $diff_in_ref_bp) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{"08.mismatch_in_5prime"}++;
                next;
            }
            substr($seq2, -$diff_in_ref_bp, $diff_in_ref_bp, "");
            print "2B $seq1\n2B $seq2\n" if ($debug);
        }


        ## ------------------------------------------------------------------------------
        ## Check that both sequences are identical on the overlap
        ## ------------------------------------------------------------------------------
        if ($seq1 ne $seq2) {
            $stats->{"09.mismatch_between_reads"}++;
            print "MM $seq1\nMM $seq2\n" if ($debug);
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Classify the PE reads into DEL, INS, WT or other-mismatches. Stores the DEL and INS in
        ## the $data_file
        ## ------------------------------------------------------------------------------
        if ($md1 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $md2 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $cigar1 =~ /^\d+M\d*D\d+M$/ and $cigar2 =~ /^\d+M\d*D\d+M$/) {
            $stats->{"OK:clean_deletion"}++;
            my ($deletion_length) = $cigar1 =~ /(\d+)D/;
#             my ($insertion_length2) = $cigar2 =~ /(\d+)D/;
#             my ($insertion_length3) = length(($md1 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#             my ($insertion_length4) = length(($md2 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#             print join("--", $insertion_length, $insertion_length2, $insertion_length3, $insertion_length4, $cigar1, $cigar2, $md1, $md2), "\n";
#             die if ($insertion_length != $insertion_length2);
#             die if ($insertion_length != $insertion_length3);
#             die if ($insertion_length != $insertion_length4);
            my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
            my $deletion_seq = substr($original_seq1, $position, $deletion_length);
            print DATA join("\t", $qname1, "DEL", $deletion_length, $position, $position + $deletion_length, $position + $deletion_length/2, $deletion_seq), "\n";
        } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M\d*I\d+M$/ and $cigar2 =~ /^\d+M\d*I\d+M$/) { 
            $stats->{"OK:clean_insertion"}++;
            my ($insertion_length) = $cigar1 =~ /(\d+)I/;
            my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
            my $insertion_seq = substr($original_seq1, $position, $insertion_length);
            print DATA join("\t", $qname1, "INS", $insertion_length, $position + 1, $position, $position + 1/2, $insertion_seq), "\n";
        } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M$/ and $cigar2 =~ /^\d+M$/) {
            $stats->{"OK:WILD-TYPE"}++;
            next;
        } else {
#             print join("\t", $cigar1, $md1, $cigar2, $md2), "\n";
            $stats->{"10.other_mismatches"}++;
            next;
        }

    }
    close(SAM);
    close(DATA);

    return ($data_file, $stats);
}


=head2 write_R_script

  Arg[1]        : string $label
  Arg[2]        : string $data_filename
  Arg[3]        : arrayref $top_sequences
  Arg[4]        : num $wt
  Arg[5]        : string $output_pdf_file
  Example       : my $R_file = write_R_script("test", "file.bam.data.txt", $top_sequences, 3123,
                  "test.pdf")
  Description   : Creates an R script to plot all the insertions and deletions in PDF, PNG and SVG
                  format
  Returns       : string with the filename with the resulting R code
  Exceptions    : 

=cut

sub write_R_script {
    my ($label, $data_file, $top_sequences, $wt, $output_pdf_file) = @_;

    my $num_lines = @$top_sequences;
    my $R_script = $data_file;
    $R_script =~ s/.data.txt$/.R/;

    open(R, ">$R_script") or die;
    print R "
data <- read.table('$data_file', header=T, row.names=1)

# =============================================================================
#  Subset data frame
# =============================================================================
#  Get new dataframes, one for deletions and one for insertions. The read name
#  was read as the row name. In this case, we skip the first column which
#  contains the type of event (either DEL or INS)
# =============================================================================
data.del <- subset(data, event=='DEL', select=2:ncol(data))
data.ins <- subset(data, event=='INS', select=2:ncol(data))


# =============================================================================
#  Number of events
# =============================================================================
#  Number of WT events is taken from the \$wt Perl variable, while the number of deletions
#  and insertions is taken from the number of rows for the deletions and insertion data frame
#  subsets respecively.
#
#  Also obtain the corresponding percentages (w.r.t. OK sequences only)
# =============================================================================
num.wt = $wt
num.del = dim(data.del)[1]
num.ins = dim(data.ins)[1]
num.total = num.wt + num.del + num.ins

perc.wt = paste0(format(100*num.wt/num.total, digits=3),'%')
perc.del = paste0(format(100*num.del/num.total, digits=3),'%')
perc.ins = paste0(format(100*num.ins/num.total, digits=3),'%')


# =============================================================================
#  FUNCTION plot.size.histograms
# =============================================================================
#  This method plots 3 histograms, one for deletions, one for insertions and
#  one with both deletions and insertions as a 'mirror' plot.
#
#  The histograms are for the lengths (sizes) of the deletions and insertions.
# =============================================================================
plot.size.histograms <- function(data.del, data.ins) {
    ### Use same breaks (x-axis) for all 3 plots. Make sure the range goes from 1-10 at least.
    breaks = (min(data.del[,1],data.ins[,1],1)-1):max(data.del[,1], data.ins[,1], 10)

    ### Get the histograms for deletions and insertions, using the set breaks. Store the
    ### result instead of plotting it.
    h.del = hist(data.del[,1], breaks=breaks, plot=F);
    h.ins = hist(data.ins[,1], breaks=breaks, plot=F);

    ### Cosmetic variables
    col.del = rgb(0.8,0,0)
    col.ins = rgb(1,0.5,0.5)
    ylim.max = max(h.del\$counts, h.ins\$counts)
    xlim = c(min(h.del\$breaks)+1,max(h.del\$breaks))
    xlab.del = paste0('Deletion sizes (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Insertion sizes (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    ### Plot histogram for deletions only (if any)
    main=paste0('Histogram of Deletion sizes (', '$label', ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    ### Plot histogram for insertions only (if any)
    main=paste0('Histogram of Insertion sizes (', '$label', ')')
    if (num.ins > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.ins\$mids, 0, h.ins\$mids+1, h.ins\$counts, col=col.ins)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }

    ### Plot histograms for both Deletion and Insertions (as a mirror plot)
    # Change the margins to make room for the upper axis + label + title
    mar = par('mar')
    par('mar' = c(mar[1]-0.5, mar[2], mar[3]+3, mar[4]))
    # Plot the data
    plot(h.del\$counts, xlim=xlim, type='n', xlab=NA, ylab=ylab,
        main=NA, ylim=c(-ylim.max, ylim.max));
    rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    rect(h.ins\$mids, 0, h.ins\$mids+1, -h.ins\$counts, col=col.ins)
    # Add the top axis
    axis(3)
    # Write the main title (as an mtext so we can have it a little higher than by default)
    mtext(paste0('Histogram of event sizes (', '$label', ')'), side=3, line = 5, cex=1.2, font=2)
    # Write the labels for the top and bottom x-axes
    mtext(xlab.del, side=3, line = 2.5)
    mtext(xlab.ins, side=1, line = 2.5)
    # Reset margins to previous value
    par('mar' = mar)
}


# =============================================================================
#  FUNCTION plot.location.images
# =============================================================================
#  This method plots 3 histograms, one for deletions, one for insertions and
#  one with both deletions and insertions as a 'mirror' plot. It then plots
#  2 scatter plots (location vs length/size of the event), one for deletions
#  and the other for insertions.
#
#  The histograms are for the location of the mid-point of the deletions and
#  insertions.
# =============================================================================
plot.location.images <- function(data.del, data.ins) {
    ### Use same breaks (x-axis) for all 3 plots. Make sure the range spans at least 10 bp.
    breaks = (min(data.del[,4],data.ins[,4])-1):max(data.del[,4], data.ins[,4])+1
    while (length(breaks) < 10) {
        if (breaks[1] > 0) {
            breaks = c(breaks[1]-1, breaks)
        }
        breaks = c(breaks, breaks[length(breaks)]+1)
    }

    ### Get the histograms for deletions and insertions, using the set breaks. Store the
    ### result instead of plotting it.
    h.del = hist(data.del[,4], breaks=breaks, plot=F);
    h.ins = hist(data.ins[,4], breaks=breaks, plot=F);

    ### Cosmetic variables
    col.del = rgb(0,0,0.8)
    col.ins = rgb(0.5,0.5,1)
    ylim.max = 1.04*max(h.del\$counts, h.ins\$counts)
    xlim = c(breaks[1],breaks[length(breaks)])
    xlab.del = paste0('Location of Deletions (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Location of Insertions (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    ### Plot histogram for deletions only (if any)
    main=paste0('Midpoint location of the deletions (', '$label', ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    ### Plot histogram for insertions only (if any)
    main=paste0('Midpoint location of the insertions (', '$label', ')')
    if (num.ins) {
        plot(h.ins\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.ins\$mids, 0, h.ins\$mids+1, h.ins\$counts, col=col.ins)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }

    ### Plot histograms for both Deletion and Insertions (as a mirror plot)
    # Change the margins to make room for the upper axis + label + title
    mar = par('mar')
    par('mar' = c(mar[1]-0.5, mar[2], mar[3]+3, mar[4]))
    # Plot the data
    plot(h.del\$counts, xlim=xlim, type='n', xlab=NA, ylab=ylab,
        main=NA, ylim=c(-ylim.max, ylim.max));
    rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    rect(h.ins\$mids, 0, h.ins\$mids+1, -h.ins\$counts, col=col.ins)
    # Add the top axis
    axis(3)
    # Write the main title (as an mtext so we can have it a little higher than by default)
    mtext(paste0('Midpoint location (', '$label', ')'), side=3, line = 5, cex=1.2, font=2)
    # Write the labels for the top and bottom x-axes
    mtext(xlab.del, side=3, line = 2.5)
    mtext(xlab.ins, side=1, line = 2.5)
    # Reset margins to previous value
    par('mar' = mar)


    ### Plot the scatter plot for deletions (if any). Uses the same x-axis as the other plots
    main = paste0('Scatter plot of Deletions sizes vs Midpoint location (', '$label', ')')
    if (num.del > 0) {
        ylim = c(0,max(data.del[,1],10))
        smoothScatter(data.del[,4], data.del[,1], xlim=xlim, ylim=ylim,
            main=main, xlab=xlab.del, ylab='Deletion size');
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }
    
    ### Plot the scatter plot for insertions (if any). Uses the same x-axis as the other plots
    main = paste0('Scatter plot of Insertion sizes vs Midpoint location (', '$label', ')')
    if (num.ins > 0) {
        ylim = c(0,max(data.ins[,1],10))
        smoothScatter(data.ins[,4], data.ins[,1], xlim=xlim, ylim=ylim,
            main=main, xlab=xlab.ins, ylab='Insertion size');
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }
}


# =============================================================================
#  FUNCTION plot.deletions.frequencies
# =============================================================================
#  This method plots a histogram showing the number of times a given basepair
#  from the wild-type sequence has been deleted. That number is a sum across
#  sequences with a clean deletion.
# =============================================================================
plot.deletion.frequencies <- function(data) {
    ### Create a new data.frame with the from and to values
    range = data.frame(from=data[,2], to=data[,3])

    ### Create a vector to store the number of times each bp has been deleted
    del <- vector(mode='numeric', length=max(data[,3]))
    for (i in 1:dim(range)[1]) { for (p in range[i,1]:range[i,2]) { del[p] <- del[p]+1 } }

    ### Genereate the plot
    # empty plot
    plot(del, xlim=c(min(range\$from-1),max(range\$to)+1), type='n',
        main=paste0('Frequency of deletion per bp (', '$label', ')'),
        xlab='Location of deleted bp', ylab='counts')
    # add light blue rectangles with white background
    rect(min(range\$from):max(range\$to)-0.5, 0, min(range\$from):max(range\$to)+0.5,
        del[min(range\$from):max(range\$to)], col='lightblue', border='white')
    # add a black line around the profile of deletion frequencies
    lines(min(range\$from):max(range\$to)-0.5,
        del[min(range\$from):max(range\$to)], type='s', col='black')
}


# =============================================================================
#  FUNCTION plot.pie.chart
# =============================================================================
#  This method plots a pie chart with the number of WT, DEL and INS sequences
# =============================================================================
plot.pie.chart <- function() {
    pie(c(num.wt, num.del, num.ins), labels=c('WT','DEL','INS'), main=paste0('Summary of events (', '$label', ')'))
    mtext(paste0('Wild-type (n = ', num.wt, '; ', perc.wt, ')'), side=1, line=0)
    mtext(paste0('Deletions (n = ', num.del, '; ', perc.del, ')'), side=1, line=1)
    mtext(paste0('Insertions (n = ', num.ins, '; ', perc.ins, ')'), side=1, line=2)
}


# =============================================================================
#  FUNCTION plot.top.seqs
# =============================================================================
#  This method plots the most common INS and DEL. The strings are extracted
#  from the \@top_sequences Perl variable
# =============================================================================
plot.topseqs <- function(data) {
    ### Use smaller top and bottom margins
    mar = par('mar')
    par('mar' = c(mar[1],1,mar[3],1))
    ### Initiate plot
    plot(NA,xlim=c(0,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA,
        main=paste0('Most common events (', '$label', ')'))
    ### Add text
    text(0, 1-1:$num_lines/25, cex=0.4, adj=0, family='mono', labels=c('", join("', '", @$top_sequences), "'))
}


# =============================================================================
#  FUNCTION plot.figures
# =============================================================================
#  This method calls all the previous methods to draw all the figures in order.
#  The advantage of having this method is that we can now easily create a PDF
#  as well as a set of PNG and SVG files (see below)
# =============================================================================
plot.figures <- function() {

    plot.pie.chart()

    ### Check that there are data to be plotted
    if (num.del+num.ins > 0) {
    
        plot.size.histograms(data.del, data.ins)
        plot.location.images(data.del, data.ins)
        plot.deletion.frequencies(data.del)

    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main='$label')
        text(0, 0, labels=c('No events to plot'))
    }
";
    if (@$top_sequences) {
        print R "   plot.topseqs(data.del)\n";
    }
    print R "
}

# =============================================================================
#  Create and save the plots
# =============================================================================
#  Here the plots are created and saved to a PDF file , then to a set of PNG
#  files and finally to a set of SVG files.
# =============================================================================
pdf('$output_pdf_file')
plot.figures();
dev.off()

# Substitute the .pdf extension by %02d.png (to create several files like 01, 02, etc)
png(sub('.pdf', '.%02d.png', '$output_pdf_file'))
plot.figures();
dev.off()

# Substitute the .pdf extension by %02d.svg (to create several files like 01, 02, etc)
svg(sub('.pdf', '.%02d.svg', '$output_pdf_file'))
plot.figures();
dev.off()

";

    close(R);
    
    return($R_script);
}
