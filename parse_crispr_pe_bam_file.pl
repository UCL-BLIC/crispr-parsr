#! /usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;

my $help;
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

Options:
--help: shows this help
--debug: output debugging information
--ref_seq: fasta reference sequence file

};

GetOptions(
    "help"  => \$help,
    "debug"  => \$debug,
    "input_file=s" => \$input_bam_file,
    "output_file=s" => \$output_R_file,
    "ref_seq=s" => \$ref_seq_file,
    );

if ($help or !$input_bam_file or !$output_R_file) {
    print $desc;
    exit();
}

open(IN_SAM, "samtools view -H $input_bam_file |") or die;
my @header = <IN_SAM>;
close(IN_SAM);

open(IN_SAM, "samtools view $input_bam_file |") or die;

my $data_file = "$input_bam_file.data.txt";
open(DATA, ">$data_file");
print DATA join("\t", "read", "deletion_length", "from", "to", "midpoint"), "\n";

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
        print DATA join("\t", $qname1, $deletion_length, $position, $position + $deletion_length, $position + $deletion_length/2), "\n";
    } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M\d*I\d+M$/ and $cigar2 =~ /^\d+M\d*I\d+M$/) { 
        $exceptions->{"OK:clean_insertion"}++;
        my ($insertion_length) = $cigar1 =~ /(\d+)I/;
        my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
        print DATA join("\t", $qname1, -$insertion_length, $position, $position + 1, $position + 1/2), "\n";
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

my $wt = $exceptions->{"OK:WILD-TYPE"};

if ($ref_seq_file) {
    open(REF, $ref_seq_file) or die;
    my $ref_seq;
    while (<REF>) {
        chomp;
        next if (/^>/);
        $ref_seq .= $_;
    }

    my @lines = qx"more $data_file | grep -v 'deletion_length' | cut -f 2,3 | grep -v -e '-' | sort | uniq -c | sort -rn | head -n 5";
    my $min_from;
    my $max_to;
    foreach my $this_line (@lines) {
        my ($num, $del_length, $from) = $this_line =~ /(\d+)\s(\d+)\s(\d+)/;
        $min_from = $from if (!$min_from or $from < $min_from);
        $max_to = $from+$del_length if (!$max_to or $from+$del_length > $max_to);
    }
    $min_from -= 21;
    $max_to += 19;
    print "\n";
    print substr($ref_seq, $min_from, $max_to-$min_from), "\t$wt\t0\tWILD-TYPE\n";
    foreach my $this_line (@lines) {
        my ($num, $del_length, $from) = $this_line =~ /(\d+)\s(\d+)\s(\d+)/;
        print substr($ref_seq, $min_from, $from-$min_from-1), '-'x$del_length, substr($ref_seq, ($from+$del_length-1), $max_to - ($from+$del_length-1)), "\t", join("\t", $num, $del_length, $from), "\n";
    }
}


open(R, "|R --vanilla --slave") or die;

print R "
data <- read.table('$data_file', header=T, row.names=1)

data <- subset(data, deletion_length > 0)

del = dim(data)[1]
wt = $wt

pdf('$output_R_file')
sub = paste0('$input_bam_file (', del, ' deletions / $wt wild-type = ', format(100*del/(del+wt), digits=3),'%)')

h = hist(data[,1], breaks=(min(data[,1])-1):max(data[,1]), plot=F);
plot(h\$counts, xlim=c(min(h\$breaks)+1,max(h\$breaks)), type='n', main='Histogram of Deletion sizes', sub=sub, xlab='Deletion size', ylab='counts');
rect(h\$mids, 0, h\$mids+1, h\$counts, col='red')

h = hist(data[,4], breaks=(min(data[,4])-1):(max(data[,4])+1), plot=F);
plot(h\$counts, xlim=c(min(h\$breaks)+1,max(h\$breaks)), type='n', main='Midpoint location of the deletion', sub=sub, xlab='Location', ylab='counts');
rect(h\$mids, 0, h\$mids+1, h\$counts, col='blue')

range = data.frame(from=data[,2], to=data[,3])
del <- vector(mode='numeric', length=max(data[,3]))
for (i in 1:dim(range)[1]) { for (p in range[i,1]:range[i,2]) { del[p] <- del[p]+1 } }
plot(del, xlim=c(min(range\$from-1),max(range\$to)+1), type='n', main='Frequency of deletion per bp', sub=sub, xlab='Location', ylab='counts')
rect(min(range\$from):max(range\$to)-0.5, 0, min(range\$from):max(range\$to)+0.5, del[min(range\$from):max(range\$to)], col='lightblue', border='white')
lines(min(range\$from):max(range\$to)-0.5, del[min(range\$from):max(range\$to)], type='s', col='black')

smoothScatter(data[,4], data[,1], main='Scatter plot of Deletions sizes vs Midpoint location' , sub=sub, xlab='Location', ylab='Deletion size');

dev.off()
";

close(R);

=pod


data <- read.table('bam/sample2.2.bam.data.txt', header=T, row.names=1)
best.hits = data.frame()
for (h in 1:5) {
  t = table(data[,1:2])
  which(t == max(t), arr.ind=T, useNames=F)
  ind <- which(t == max(t), arr.ind=T, useNames=F)
  next.hit = c(as.numeric(rownames(t)[ind[1]]), as.numeric(colnames(t)[ind[2]]))
  next.hit = c(next.hit, from=next.hit[2]-next.hit[1]/2, to=next.hit[2]+next.hit[1]/2, num = dim(subset(data, deletion_length==next.hit[1] & location==next.hit[2]))[1])
  if (next.hit[5] != max(t)) {
    break
  }
#  print(subset(data, deletion_length==next.hit[1] & location==next.hit[2]))
#  print(dim(data))
  data <- subset(data, deletion_length!=next.hit[1] | location!=next.hit[2])
#  print(dim(data))
#  print(subset(data, deletion_length==next.hit[1] & location==next.hit[2]))
  print(next.hit)
  best.hits = rbind(best.hits, next.hit)
}
colnames(best.hits) <- c('deletion_length', 'location', 'from', 'to', 'num')
print(best.hits)
seq.from = min(best.hits$from)-20
seq.to = max(best.hits$to)+20
print(paste(seq.from, seq.to))
print(paste0(rep('N', seq.to - seq.from + 1), collapse=""))
for (i in 1:nrow(best.hits)) {
  print(paste0(
    paste0(rep('N', best.hits[i,3] - seq.from), collapse=""),
    paste0(rep('-', best.hits[i,4] - best.hits[i,3] + 1), collapse=""),
    paste0(rep('N', seq.to - best.hits[i,4]), collapse=""),
    '  ', best.hits[i,1], '  ', best.hits[i,3], '  ', best.hits[i,4], collapse=""))
}



[1] "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
[1] "NNNNNNNNNNNNNNNNNNNNNNN------------------NNNNNNNNNNNNNNNNNNNNNNNNN  17  118  135"
[1] "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-----NNNNNNNNNNNNNNNNNNNNNNNNNNNNN  4  127  131"
[1] "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN-------NNNNNNNNNNNNNNNNNNNNNNNNNNNNN  6  125  131"
[1] "NNNNNNNNNNNNNNNNNNNN--------------------------NNNNNNNNNNNNNNNNNNNN  25  115  140"
[1] "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN----NNNNNNNNNNNNNNNNNNNNNNNNNNNNN  3  128  131"



ATGATGTTTTCCCCCCGTTTTTGTTTTTTGTTTTGTAGTTGATATTCACTGATGGACTCCTGGTAGAGAAGAAAACCCCAGCAGTGTGCTTGCTCAGGAGAGGGGAGATGTGATGGA
TTGATATTCACTGATGGACTCCAAAGAATCATTAACTCCTGGTAGAGAAGAAAACCCCAGCAGTGTGCTTGCTCAGGAGA
TTGATATTCACTGATGGAACTCCAAAGAATCATTACTCCTGGTAGAGAAGAAAACCCCAGCAGTGTGCTTGCTCAGGAGAGGGGAGATGTGATGGACTTCTATAAAACCCTAAGAGGAGGAGCTACTGTG
ttgaTATTCACTGATGGA-----------------CTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc
ttgaTATTCACTGATGG-----------------CTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc
ttgaTATTCACTGATGGACTCCAAAGAATCATTACTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc
TTGATATTCACTGATGGACTCCAAAGAATCATTAACTCCTGGTAGAGAAGAAAACCCCAGCAGTGTGCTTGCTCAGGAGA




ttgaTATTCACTGATGGACTCCAAAGAATCATTAACTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc
ttgaTATTCACTGATGG-----------------CTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc	472	17	118
ttgaTATTCACTGATGGACTCCAAAG----TTAACTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc	110	4	127
ttgaTATTCACTGATGGACTCCAA------TTAACTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc	55	6	125
ttgaTATTCACTGA-------------------------GGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc	52	25	115
ttgaTATTCACTGATGGACTCCAAAGA---TTAACTCCTGGTAGAGAAGAAAaccccagcagtgtgcttgctcaggagaggggagatgtgatggacttctataaaaccc	50	3	128
