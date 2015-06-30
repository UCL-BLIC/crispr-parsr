#! /usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $help;
my $label;
my $input_bam_file;
my $output_R_file;
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
    "output_file=s" => \$output_R_file,
    "ref_seq=s" => \$ref_seq_file,
    );

if ($help or !$input_bam_file or !$output_R_file) {
    print $desc;
    exit();
}

if (!$label) {
    $label = $input_bam_file;
}

open(IN_SAM, "samtools view -H $input_bam_file |") or die;
my @header = <IN_SAM>;
close(IN_SAM);

open(IN_SAM, "samtools view $input_bam_file |") or die;

my $data_file = "$input_bam_file.data.txt";
open(DATA, ">$data_file");
print DATA join("\t", "read", "event", "event_length", "from", "to", "midpoint", "seq"), "\n";

my $exceptions;
while (<IN_SAM>) {
    chomp;
    my ($qname1, $flag1, $rname1, $pos1, $mapq1, $cigar1, $rnext1, $pnext1, $tlen1, $seq1, $qual1, @others1) = split("\t", $_);
    $_ = <IN_SAM>;
    chomp;
    my ($qname2, $flag2, $rname2, $pos2, $mapq2, $cigar2, $rnext2, $pnext2, $tlen2, $seq2, $qual2, @others2) = split("\t", $_);
    my $md1 = (grep {/^MD:Z:/} @others1)[0];
    my $md2 = (grep {/^MD:Z:/} @others2)[0];
    if ($debug) {
      print join("\t", $qname1, $flag1, $rname1, $pos1, $mapq1, $cigar1, $rnext1, $pnext1, $tlen1), "\n";
      print join("\t", $qname2, $flag2, $rname2, $pos2, $mapq2, $cigar2, $rnext2, $pnext2, $tlen2), "\n";
    }
    if ($qname1 ne $qname2) {
        die "SAM file sorted or not for PE reads\n";
    }
    next if ($rname1 ne $rname2);
    if ($cigar1 eq "*" or $cigar2 eq "*") {
        $exceptions->{"01.no_cigar"}++;
        next;
    }

    if (length($seq1) < $min_length or length($seq2) < $min_length) {
        $exceptions->{"02.short_read"}++;
        next;
    }

    my $end1 = $pos1 + length($seq1);
    my $end2 = $pos2 + length($seq2);
    my $overlap_start = $pos1>$pos2?$pos1:$pos2;
    my $overlap_end = $end1<$end2?$end1:$end2;

    if ($overlap_end < $overlap_start + $min_overlap) {
        $exceptions->{"03.no_overlap"}++;
        next;
    }
#     print "Overlap: $overlap_start-$overlap_end\n";

    my $original_seq1 = $seq1;
    my $original_seq2 = $seq2;

    if ($pos1 > $pos2) {
        # Clip R2
        my $diff_in_ref_bp = $pos1 - $pos2;
        my ($initial_match) = $cigar2 =~ /^(\d*)M/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
        if ($initial_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"05.indel_in_3prime"}++;
            next;
        }
        ($initial_match) = $md2 =~ /^MD:Z:(\d*)/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
        if ($initial_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"06.mismatch_in_3prime"}++;
            next;
        }
        substr($seq2, 0, $diff_in_ref_bp, "");
        print "1A $seq1\n1A $seq2\n" if ($debug);
    } elsif ($pos1 < $pos2) {
        # Clip R1
        my $diff_in_ref_bp = $pos2 - $pos1;
        my ($initial_match) = $cigar1 =~ /^(\d*)M/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
        if ($initial_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"05.indel_in_3prime"}++;
            next;
        }
        ($initial_match) = $md1 =~ /^MD:Z:(\d*)/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
        if ($initial_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"06.mismatch_in_3prime"}++;
            next;
        }
        substr($seq1, 0, $diff_in_ref_bp, "");
        print "1B $seq1\n1B $seq2\n" if ($debug);
    }

    if ($end1 > $end2) {
        # Clip R1
        my $diff_in_ref_bp = $end1 - $end2;
        my ($last_match) = $cigar1 =~ /(\d*)M$/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
        if ($last_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"07.indel_in_5prime"}++;
            next;
        }
        ($last_match) = $md1 =~ /(\d*)$/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
        if ($last_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"08.mismatch_in_5prime"}++;
            next;
        }
        substr($seq1, -$diff_in_ref_bp, $diff_in_ref_bp, "");
        print "2A $seq1\n2A $seq2\n" if ($debug);
    } elsif ($end1 < $end2) {
        # Clip R2
        my $diff_in_ref_bp = $end2 - $end1;
        my ($last_match) = $cigar2 =~ /(\d*)M$/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
        if ($last_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"07.indel_in_5prime"}++;
            next;
        }
        ($last_match) = $md2 =~ /(\d*)$/;
#         print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
        if ($last_match < $diff_in_ref_bp) {
            # Indels in the 3' non-overlapping sequence: ignore the pair!
            $exceptions->{"08.mismatch_in_5prime"}++;
            next;
        }
        substr($seq2, -$diff_in_ref_bp, $diff_in_ref_bp, "");
        print "2B $seq1\n2B $seq2\n" if ($debug);
    }

    if ($seq1 ne $seq2) {
        $exceptions->{"09.mismatch_in_overlap"}++;
        print "MM $seq1\nMM $seq2\n" if ($debug);
        next;
    }

    if ($md1 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $md2 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $cigar1 =~ /^\d+M\d*D\d+M$/ and $cigar2 =~ /^\d+M\d*D\d+M$/) {
        $exceptions->{"OK:clean_deletion"}++;
        my ($deletion_length) = $cigar1 =~ /(\d+)D/;
#         my ($insertion_length2) = $cigar2 =~ /(\d+)D/;
#         my ($insertion_length3) = length(($md1 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#         my ($insertion_length4) = length(($md2 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#         print join("--", $insertion_length, $insertion_length2, $insertion_length3, $insertion_length4, $cigar1, $cigar2, $md1, $md2), "\n";
#         die if ($insertion_length != $insertion_length2);
#         die if ($insertion_length != $insertion_length3);
#         die if ($insertion_length != $insertion_length4);
        my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
        my $deletion_seq = substr($original_seq1, $position, $deletion_length);
        print DATA join("\t", $qname1, "DEL", $deletion_length, $position, $position + $deletion_length, $position + $deletion_length/2, $deletion_seq), "\n";
    } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M\d*I\d+M$/ and $cigar2 =~ /^\d+M\d*I\d+M$/) { 
        $exceptions->{"OK:clean_insertion"}++;
        my ($insertion_length) = $cigar1 =~ /(\d+)I/;
        my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
        my $insertion_seq = substr($original_seq1, $position, $insertion_length);
        print DATA join("\t", $qname1, "INS", $insertion_length, $position + 1, $position, $position + 1/2, $insertion_seq), "\n";
    } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M$/ and $cigar2 =~ /^\d+M$/) {
        $exceptions->{"OK:WILD-TYPE"}++;
        next;
    } else {
#         print join("\t", $cigar1, $md1, $cigar2, $md2), "\n";
        $exceptions->{"10.other_mismatches"}++;
        next;
    }
}
close(IN_SAM);
close(DATA);

my $total = 0;
foreach my $key (sort keys %$exceptions) {
    my $value = $exceptions->{$key};
    print "$key: $value\n";
    $total += $value;
}
print "TOTAL: $total\n";

my $wt = ($exceptions->{"OK:WILD-TYPE"} or 0);

my $wt_sequence;
my @top_sequences;
if ($ref_seq_file) {
    open(REF, $ref_seq_file) or die;
    my $ref_seq;
    while (<REF>) {
        chomp;
        next if (/^>/);
        $ref_seq .= $_;
    }

    my @del_lines = qx"more $data_file | awk '\$2 == \"DEL\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn | head -n 10";
    my @ins_lines = qx"more $data_file | awk '\$2 == \"INS\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn | head -n 10";
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
    $max_to += 19;

    my $format = "\%-".($max_to-$min_from)."s \%7s \%4s \%3s \%6s \%-${longest_seq}s";
    my $header = sprintf($format, "Sequence" , "Num", "TYPE", "L", "POS", "Diff");
    my $wt_sequence = sprintf($format, substr($ref_seq, $min_from, $max_to-$min_from), $wt, "W-T", 0, "NA", "");

    @top_sequences = ($header, "", $wt_sequence, "");

    foreach my $this_line (@del_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $from-$min_from-1) . '-'x$del_length . substr($ref_seq, ($from+$del_length-1), $max_to - ($from+$del_length-1));
        push(@top_sequences, sprintf($format, $aligned_seq, $num, "DEL", $del_length, $from, $seq));
    }
    push(@top_sequences, "");
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $max_to-$min_from);
        substr($aligned_seq, $from-$min_from-2, 2, "><");
        push(@top_sequences, sprintf($format, $aligned_seq, $num, "INS", $ins_length, $from, ">$seq<"));
    }
    

    print "\n", join("\n", @top_sequences), "\n";
}

my $num_lines = @top_sequences;

open(R, "|R --vanilla --slave") or die;

print R "
data <- read.table('$data_file', header=T, row.names=1)

data.del <- subset(data, event=='DEL', select=2:ncol(data))

data.ins <- subset(data, event=='INS', select=2:ncol(data))

num.wt = $wt
num.del = dim(data.del)[1]
num.ins = dim(data.ins)[1]
num.total = num.wt + num.del + num.ins

perc.wt = paste0(format(100*num.wt/num.total, digits=3),'%')
perc.del = paste0(format(100*num.del/num.total, digits=3),'%')
perc.ins = paste0(format(100*num.ins/num.total, digits=3),'%')

plot.size.histograms <- function(data.del, data.ins, sub) {
    col.del = rgb(0.8,0,0)
    col.ins = rgb(1,0.5,0.5)
    breaks = (min(data.del[,1],data.ins[,1],1)-1):max(data.del[,1], data.ins[,1], 10)

    h.del = hist(data.del[,1], breaks=breaks, plot=F);
    h.ins = hist(data.ins[,1], breaks=breaks, plot=F);
    ylim.max = max(h.del\$counts, h.ins\$counts)
    xlim = c(min(h.del\$breaks)+1,max(h.del\$breaks))
    xlab.del = paste0('Deletion sizes (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Insertion sizes (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    # Deletions only
    main=paste0('Histogram of Deletion sizes (', '$label', ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    # Insertions only
    main=paste0('Histogram of Insertion sizes (', '$label', ')')
    if (num.ins > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.ins\$mids, 0, h.ins\$mids+1, h.ins\$counts, col=col.ins)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }

    # Both Deletion and Insertions (as a mirror plot)
    mar = par('mar')
    par('mar' = c(mar[1]-0.5, mar[2], mar[3]+3, mar[4]))
    plot(h.del\$counts, xlim=xlim, type='n', xlab=NA, ylab=ylab,
        main=NA, ylim=c(-ylim.max, ylim.max));
    axis(3)
    mtext(paste0('Histogram of event sizes (', '$label', ')'), side=3, line = 5, cex=1.2, font=2)
    mtext(xlab.del, side=3, line = 2.5)
    mtext(xlab.ins, side=1, line = 2.5)
    rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    rect(h.ins\$mids, 0, h.ins\$mids+1, -h.ins\$counts, col=col.ins)
    par('mar' = mar)
}

plot.location.histograms <- function(data.del, data.ins, sub) {
    col.del = rgb(0,0,0.8)
    col.ins = rgb(0.5,0.5,1)
    breaks = (min(data.del[,4],data.ins[,4])-1):max(data.del[,4], data.ins[,4])+1
    while (length(breaks) < 10) {
        if (breaks[1] > 0) {
            breaks = c(breaks[1]-1, breaks)
        }
        breaks = c(breaks, breaks[length(breaks)]+1)
    }
    h.del = hist(data.del[,4], breaks=breaks, plot=F);
    h.ins = hist(data.ins[,4], breaks=breaks, plot=F);
    ylim.max = 1.04*max(h.del\$counts, h.ins\$counts)
    xlim = c(breaks[1],breaks[length(breaks)])
    xlab.del = paste0('Location of Deletions (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Location of Insertions (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    # Deletions only
    main=paste0('Midpoint location of the deletions (', '$label', ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    # Insertions only
    main=paste0('Midpoint location of the insertions (', '$label', ')')
    if (num.ins) {
        plot(h.ins\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.ins\$mids, 0, h.ins\$mids+1, h.ins\$counts, col=col.ins)
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }

    # Both Deletion and Insertions (as a mirror plot)
    mar = par('mar')
    par('mar' = c(mar[1]-0.5, mar[2], mar[3]+3, mar[4]))
    plot(h.del\$counts, xlim=xlim, type='n', xlab=NA, ylab=ylab,
        main=NA, ylim=c(-ylim.max, ylim.max));
    axis(3)
    mtext(paste0('Midpoint location (', '$label', ')'), side=3, line = 5, cex=1.2, font=2)
    mtext(xlab.del, side=3, line = 2.5)
    mtext(xlab.ins, side=1, line = 2.5)
    rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    rect(h.ins\$mids, 0, h.ins\$mids+1, -h.ins\$counts, col=col.ins)
    par('mar' = mar)

    main = paste0('Scatter plot of Deletions sizes vs Midpoint location (', '$label', ')')
    if (num.del > 0) {
        smoothScatter(data.del[,4], data.del[,1], xlim=xlim,
            main=main, xlab=xlab.del, ylab='Deletion size');
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }
    
    main = paste0('Scatter plot of Insertion sizes vs Midpoint location (', '$label', ')')
    if (num.ins > 0) {
        smoothScatter(data.ins[,4], data.ins[,1], xlim=xlim,
            main=main, xlab=xlab.ins, ylab='Insertion size');
    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }
}

plot.deletion.frequencies <- function(data, sub) {
    range = data.frame(from=data[,2], to=data[,3])
    del <- vector(mode='numeric', length=max(data[,3]))
    for (i in 1:dim(range)[1]) { for (p in range[i,1]:range[i,2]) { del[p] <- del[p]+1 } }
    plot(del, xlim=c(min(range\$from-1),max(range\$to)+1), type='n',
        main=paste0('Frequency of deletion per bp (', '$label', ')'),
        sub=sub, xlab='Location', ylab='counts')
    rect(min(range\$from):max(range\$to)-0.5, 0, min(range\$from):max(range\$to)+0.5,
        del[min(range\$from):max(range\$to)], col='lightblue', border='white')
    lines(min(range\$from):max(range\$to)-0.5,
        del[min(range\$from):max(range\$to)], type='s', col='black')
}

plot.pie.chart <- function() {
    pie(c(num.wt, num.del, num.ins), labels=c('WT','DEL','INS'), main=paste0('Summary (', '$label', ')'))
    mtext(paste0('Wild-type (n = ', num.wt, '; ', perc.wt, ')'), side=1, line=0)
    mtext(paste0('Deletions (n = ', num.del, '; ', perc.del, ')'), side=1, line=1)
    mtext(paste0('Insertions (n = ', num.ins, '; ', perc.ins, ')'), side=1, line=2)
}

plot.topseqs <- function(data, sub) {
    mar = par('mar')
    par('mar' = c(mar[1],1,mar[3],1))
    plot(NA,xlim=c(0,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA,
        main=paste0('Most common events (', '$label', ')'))
    text(0, 1-1:$num_lines/25, cex=0.4, adj=0, family='mono', labels=c('", join("', '", @top_sequences), "'))
}

plot.figures <- function() {
    if (num.del+num.ins > 0) {
    
        sub = paste0(num.del, ' deletions / $wt wild-type = ', format(100*num.del/(num.del+num.wt), digits=3),'%')

        plot.size.histograms(data.del, data.ins, sub)
        plot.location.histograms(data.del, data.ins, sub)
        plot.deletion.frequencies(data.del, sub)
        plot.pie.chart()

    } else {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main='$label')
        text(0, 0, labels=c('No events to plot'))
    }
    plot.topseqs(data.del, sub)

}

pdf('$output_R_file')
plot.figures();
dev.off()

png(sub('.pdf', '.%02d.png', '$output_R_file'))
plot.figures();
dev.off()

svg(sub('.pdf', '.%02d.svg', '$output_R_file'))
plot.figures();
dev.off()

";

close(R);


