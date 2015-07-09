#! /usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;

=pod

=head1 NAME

crispr-plotr - Software to plot deletions and insertions from a CRISPR resequencing experiment

=head1 SYNOPSIS

crispr-plotr.pl [--ref_seq ref_seq.fa] [--label Sample31] --input crispr.bam --output Sampel31.pdf

=head1 DESCRIPTION

Briefly, this software takes a bam file with overlapping paired-end (PE) reads and call insertions
and deletions with respect to the reference sequence on the overlap.

From the input set of PE reads, it will filter out any pair where:

=over 8

=item * any of the reads does not map

=item * any of the reads is too short

=item * the reads do not overlap significantly

=item * have mutations or indels before or after the overlap

=item * have more than just a deletion or just an insertion in the overlap

=back

The output is a series of plots with overall stats on the number of PE reads that pass or not the
filters; the proportion of WT, insertions and deletions; the list of most common events; the
distributions of insertion and deletion sizes as well as the distribution of the locations of these
events; a scatter plot showing the relationship between both and one last plot with the frequency
of deletions per bp.

These plots are generated with R. The R script required to generate the plots is also stored in
the output folder and can be edited by the user to modify the plots if needed.

=head1 OPTIONS

Please note that options names may be abbreviated to uniqueness, case does not matter, and a
single dash is sufficient

=over 8

=item B<-h>|B<--help>

Prints the full help.

=item B<-i>|B<--input crispr.bam>

The BAM filename with the CRISPR aligned PE reads.

=item B<-o>|B<--output FILE.pdf>

By default, the output is PDF file. Other output are available and the filename will the same as
this one, but with the required extension.

=item B<--label LABEL>

Use this as the sample name for labelling the plots.

=item B<--min-overlap 10>

Minimum overlap between the PE reads. Pairs that do not overlap by at least this length will be
filtered out.

=item B<--min-length 80>

Minimum length of the reads. If either read is shorter, the pair is filtered out

=item B<--allow-indels>

Do not filter out PE read that have indels in their flanking region, where the flanking region is
the region outside of the overlap, i.e. covered by one of the reads only.

=item B<--allow-muts>

Do not filter out PE read that have mutations in their flanking region, where the flanking region is
the region outside of the overlap, i.e. covered by one of the reads only.

=item B<--allow-any>

Sets the B<allow-indels> and B<allow-muts> flags to true.

=item B<--ref-seq>|B<--wt-seq> ref_seq.fa>

This is a simple FASTA file with one sequence only (typically a short one) corresponding to the expected
wild-type sequence, before editing has ocurred.

    Example:
    ---------------------------------------
    >sample_amplicon1_ref_123456
    AACAGTGGGTCTTCGGGAACCAAGCAGCAGCCATGGGAGGTCTGGCTGTGTTCAGGCT
    CTGCTCGTGTAGATTCACAGCGCGCTCTGAACCCCCGCTGAGCTACCGATGGAAGAGG
    AGGAGGTCCTACAGTCGGAGATTCACAGCGCGCTCTGAACCACTTTCAGGAGACTCGA
    CTATTATGACTTATACGCGATA
    ---------------------------------------

If not provided, the list of most common events is not displayed.

=item B<--plot-pdf>|B<--pdf> or B<--noplot-pdf>|B<--nopdf>

Enable or disable the PDF output.

Default: PDF output is enabled.

=item B<--plot-png>|B<--png> or B<--plot-nopng>|B<--nopng>

Enable or disable the PNG output.

Default: PNG output is disabled.

=item B<--plot-svg>|B<--svg> or B<--noplot-svg>|B<--nosvg>

Enable or disable the SVG output.

Default: SVG output is disabled.

=item B<--force-plots>

This options controls whether you want to obtain always the same number of plots or not. By default,
crispr-plotr skips the plots for which there is no data. For instance, the plot for insertions will
not be generated if there are no insertions in the input BAM file.

By using this option, you are requesting to always get an image even if there are no data.
Typically the plot will contain a simple message like "No insertions". This is useful if you wish
to embed the images in an automatically generated report for instance.

=back

=head1 Requirements

=over 8

=item B<R>: http://www.r-project.org

=back

=head1 INTERNAL METHODS

The rest of the documentation refers to the internal methods within this software and is
intended for developers only.

=cut

my $help;
my $label;
my $input_bam_file;
my $output_pdf_file;
my $debug;
my $ref_seq_file;
my $min_overlap = 10;
my $min_length = 80;
my $allow_indels = 0;
my $allow_muts = 0;
my $allow_any = 0;
my $plot_pdf = 1;
my $plot_png = 0;
my $plot_svg = 0;
my $plot_all = 0;
my $force_plots = 0;

my $STAT_NO_ALIGNMENT = "No alignment";
my $STAT_SHORT_READ = "Short read";
my $STAT_NO_OVERLAP = "No overlap";
my $STAT_INDEL_5_PRIME = "5' indel";
my $STAT_MUTATION_5_PRIME = "5' mut";
my $STAT_INDEL_3_PRIME = "3' indel";
my $STAT_MUTATION_3_PRIME = "3' mut";
my $STAT_READ_MISMATCH = "Reads differ";
my $STAT_OTHER_MISMATCHES = "Other mismatches";
my $STAT_OK_WILD_TYPE = "WT";
my $STAT_OK_DELETION = "DEL";
my $STAT_OK_INSERTION = "INS";


my @STAT_ORDER = (
    $STAT_NO_ALIGNMENT,
    $STAT_SHORT_READ,
    $STAT_NO_OVERLAP,
    $STAT_INDEL_5_PRIME,
    $STAT_MUTATION_5_PRIME,
    $STAT_INDEL_3_PRIME,
    $STAT_MUTATION_3_PRIME,
    $STAT_READ_MISMATCH,
    $STAT_OTHER_MISMATCHES,
    $STAT_OK_WILD_TYPE,
    $STAT_OK_DELETION,
    $STAT_OK_INSERTION,
);

my $COLOR_DELETION = "cyan3";
my $COLOR_INSERTION = "darkorchid3";
my $STAT_COLOR = {
    "$STAT_NO_ALIGNMENT" => "lightgrey",
    "$STAT_SHORT_READ" => "grey",
    "$STAT_NO_OVERLAP" => "darkgrey",
    "$STAT_INDEL_5_PRIME" => "goldenrod1",
    "$STAT_MUTATION_5_PRIME" => "goldenrod2",
    "$STAT_INDEL_3_PRIME" => "goldenrod3",
    "$STAT_MUTATION_3_PRIME" => "goldenrod4",
    "$STAT_READ_MISMATCH" => "coral2",
    "$STAT_OTHER_MISMATCHES" => "coral3",
    "$STAT_OK_WILD_TYPE" => "white",
    "$STAT_OK_DELETION" => $COLOR_DELETION,
    "$STAT_OK_INSERTION" => $COLOR_INSERTION,
};


GetOptions(
    "help"  => \$help,
    "debug"  => \$debug,
    "label=s" => \$label,
    "input_file|input-file=s" => \$input_bam_file,
    "output_file|output-file=s" => \$output_pdf_file,
    "ref_seq|ref-seq|wt_file|wt-file|wt_seq|wt-seq=s" => \$ref_seq_file,
    "min_overlap|min-overlap=i" => \$min_overlap,
    "min_length|min-length=i" => \$min_length,

    "allow_indels|allow-indels|indels!" => \$allow_indels,
    "allow_muts|allow-muts|muts!" => \$allow_muts,
    "allow_any|allow-any|any!" => \$allow_any,

    "plot_pdf|plot-pdf|pdf!" => \$plot_pdf,
    "plot_png|plot-png|png!" => \$plot_png,
    "plot_svg|plot-svg|svg!" => \$plot_svg,
    "plot_all|plot-all!" => \$plot_all,
    "force_plots|force-plots!" => \$force_plots,
    );

if ($help) {
    pod2usage(-verbose=>2);
}

if (!$input_bam_file or !$output_pdf_file) {
    pod2usage(-verbose=>1);
}

if (!$label) {
    $label = $input_bam_file;
}

if ($allow_any) {
    $allow_indels = 1;
    $allow_muts = 1;
}

if ($plot_all) {
    $plot_pdf = 1;
    $plot_png = 1;
    $plot_svg = 1;
}

my ($data_file, $stats) = parse_bam_file($input_bam_file);

my $total = 0;
foreach my $key (@STAT_ORDER) {
    my $value = ($stats->{$key} or 0);
    print "$key: $value\n";
    $total += $value;
}
print "TOTAL: $total\n";

my $top_sequences = [];
my $ref_seq = "";
if ($ref_seq_file) {
    $ref_seq = read_fasta_seq($ref_seq_file);

    my $wt = ($stats->{$STAT_OK_WILD_TYPE} or 0);
    $top_sequences = get_top_sequences($ref_seq, $data_file, $wt);
}

my $R_script = write_R_script($label, $data_file, $ref_seq, $top_sequences, $stats, $output_pdf_file);

print qx"Rscript $R_script";

exit(0);


=head2 read_fasta_seq

  Arg[1]        : string $fasta_file
  Example       : my $seq = read_fasta_seq($fasta_file);
  Description   : Reads the sequence in the FASTA file $fasta_file
  Returns       : string $seq
  Exceptions    : Dies if any filename is not found.

=cut

sub read_fasta_seq {
    my ($fasta_file) = @_;
    my $seq;

    open(FASTA, $fasta_file) or die;
    while (<FASTA>) {
        chomp;
        next if (/^>/);
        $seq .= $_;
    }

    return $seq;
}

=head2 get_top_sequences

  Arg[1]        : string $ref_seq
  Arg[2]        : string $data_filename
  Arg[3]        : integer $num_of_wild_type_seqs
  Example       : my $top_sequences = get_top_sequences($ref_seq, $$data_file, 1213);
  Description   : Reads from the $data_file the most common deletions and insertions (up to 10) and
                  aligns them to the wild-type sequence.
  Returns       : arrayref of strings
  Exceptions    : 

=cut

sub get_top_sequences {
    my ($ref_seq, $data_file, $wt) = @_;
    my $top_sequences = [];


    ## ------------------------------------------------------------------------------
    ## Reads and extract stats on most common deletions and insertions:
    ## ------------------------------------------------------------------------------
    my @del_lines = qx"more $data_file | awk '\$2 == \"DEL\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn";
    my @ins_lines = qx"more $data_file | awk '\$2 == \"INS\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn";


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
    if (!defined($max_to)) {
        # i.e. no deletions nor insertions
        $max_to = int(length($ref_seq)/2)+10;
        $min_from = $max_to-10;
    }
    $min_from -= 21;
    $min_from = 0 if ($min_from < 0);
    $max_to += 19;

    ## ------------------------------------------------------------------------------
    ## Sets the output format (using whitespaces to be nicely printed in R afterewards)
    ## ------------------------------------------------------------------------------
    my $format = "\%-".($max_to-$min_from)."s \%7s  \%-4s \%3s \%6s \%-${longest_seq}s";


    ## ------------------------------------------------------------------------------
    ## Header and WT sequence
    ## ------------------------------------------------------------------------------
    my $header = sprintf($format, "Sequence" , "Num", "TYPE", "L", "POS", "Diff");
    my $wt_sequence = sprintf($format, substr($ref_seq, $min_from, $max_to-$min_from), $wt, "WT", 0, "NA", "");
    $top_sequences = [$header, "", $wt_sequence, ""];


    ## ------------------------------------------------------------------------------
    ## Most common deletions
    ## ------------------------------------------------------------------------------
    my @del_sequences;
    foreach my $this_line (@del_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $from-$min_from-1) . '-'x$del_length . substr($ref_seq, ($from+$del_length-1), $max_to - ($from+$del_length-1));
        push(@del_sequences, sprintf($format, $aligned_seq, $num, "DEL", $del_length, $from, $seq));
    }
    my $deletions_file = $data_file;
    $deletions_file =~ s/\.txt$/.del.txt/;
    open(DEL, ">$deletions_file") or die;
    print DEL join("\n", $header, $wt_sequence, "", @del_sequences, "");
    close(DEL);
    push(@$top_sequences, splice(@del_sequences, 0, 20));


    ## ------------------------------------------------------------------------------
    ## Separation line
    ## ------------------------------------------------------------------------------
    push(@$top_sequences, "");


    ## ------------------------------------------------------------------------------
    ## Most common insertions
    ## ------------------------------------------------------------------------------
    my @ins_sequences;
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from, $max_to-$min_from);
        substr($aligned_seq, $from-$min_from-2, 2, "><");
        push(@ins_sequences, sprintf($format, $aligned_seq, $num, "INS", $ins_length, $from, ">$seq<"));
    }
    my $insertions_file = $data_file;
    $insertions_file =~ s/\.txt$/.ins.txt/;
    open(INS, ">$insertions_file") or die;
    print INS join("\n", $header, $wt_sequence, "", @ins_sequences, "");
    close(INS);
    push(@$top_sequences, splice(@ins_sequences, 0, 20));

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
            $stats->{$STAT_NO_ALIGNMENT}++;
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Check that both reads are long enough
        ## ------------------------------------------------------------------------------
        if (length($seq1) < $min_length or length($seq2) < $min_length) {
            $stats->{$STAT_SHORT_READ}++;
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Check that both reads overlap and that the overlap is long enough
        ## ------------------------------------------------------------------------------
        my $overlap_start = $pos1>$pos2?$pos1:$pos2;
        my $overlap_end = $end1<$end2?$end1:$end2;
        if ($overlap_end < $overlap_start + $min_overlap) {
            $stats->{$STAT_NO_OVERLAP}++;
            next;
        }
#        print "Overlap: $overlap_start-$overlap_end\n";


        ## ------------------------------------------------------------------------------
        ## Keeps a safe copy of both original sequences
        ## ------------------------------------------------------------------------------
        my $original_seq1 = $seq1;
        my $original_seq2 = $seq2;


        ## ------------------------------------------------------------------------------
        ## Clip the sequence 5' of the overlap (and check that there is no mismatch nor indel)
        ## ------------------------------------------------------------------------------
        if ($pos1 > $pos2) {
            # Clip R2
            my $diff_in_ref_bp = $pos1 - $pos2;
            my ($initial_match) = $cigar2 =~ /^(\d*)M/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if (!$allow_indels and $initial_match < $diff_in_ref_bp) {
                # Indels in the 5' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_INDEL_5_PRIME}++;
                next;
            }
            ($initial_match) = $md2 =~ /^MD:Z:(\d*)/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if (!$allow_muts and $initial_match < $diff_in_ref_bp) {
                # Indels in the 5' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_MUTATION_5_PRIME}++;
                next;
            }
            substr($seq2, 0, $diff_in_ref_bp, "");
            print "1A $seq1\n1A $seq2\n" if ($debug);
        } elsif ($pos1 < $pos2) {
            # Clip R1
            my $diff_in_ref_bp = $pos2 - $pos1;
            my ($initial_match) = $cigar1 =~ /^(\d*)M/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if (!$allow_indels and $initial_match < $diff_in_ref_bp) {
                # Indels in the 5' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_INDEL_5_PRIME}++;
                next;
            }
            ($initial_match) = $md1 =~ /^MD:Z:(\d*)/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $initial_match), "\n";
            if (!$allow_muts and $initial_match < $diff_in_ref_bp) {
                # Indels in the 5' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_MUTATION_5_PRIME}++;
                next;
            }
            substr($seq1, 0, $diff_in_ref_bp, "");
            print "1B $seq1\n1B $seq2\n" if ($debug);
        }


        ## ------------------------------------------------------------------------------
        ## Clip the sequence 3' of the overlap (and check that there is no mismatch nor indel)
        ## ------------------------------------------------------------------------------
        if ($end1 > $end2) {
            # Clip R1
            my $diff_in_ref_bp = $end1 - $end2;
            my ($last_match) = $cigar1 =~ /(\d*)M$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if (!$allow_indels and ($last_match < $diff_in_ref_bp)) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_INDEL_3_PRIME}++;
                next;
            }
            ($last_match) = $md1 =~ /(\d*)$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if (!$allow_muts and ($last_match < $diff_in_ref_bp)) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_MUTATION_3_PRIME}++;
                next;
            }
            substr($seq1, -$diff_in_ref_bp, $diff_in_ref_bp, "");
            print "2A $seq1\n2A $seq2\n" if ($debug);
        } elsif ($end1 < $end2) {
            # Clip R2
            my $diff_in_ref_bp = $end2 - $end1;
            my ($last_match) = $cigar2 =~ /(\d*)M$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if (!$allow_indels and ($last_match < $diff_in_ref_bp)) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_INDEL_3_PRIME}++;
                next;
            }
            ($last_match) = $md2 =~ /(\d*)$/;
#             print join("\t", $pos1, $pos2, $diff_in_ref_bp, $last_match), "\n";
            if (!$allow_muts and ($last_match < $diff_in_ref_bp)) {
                # Indels in the 3' non-overlapping sequence: ignore the pair!
                $stats->{$STAT_MUTATION_3_PRIME}++;
                next;
            }
            substr($seq2, -$diff_in_ref_bp, $diff_in_ref_bp, "");
            print "2B $seq1\n2B $seq2\n" if ($debug);
        }


        ## ------------------------------------------------------------------------------
        ## Check that both sequences are identical on the overlap
        ## ------------------------------------------------------------------------------
        if ($seq1 ne $seq2) {
            $stats->{$STAT_READ_MISMATCH}++;
            print "MM $seq1\nMM $seq2\n" if ($debug);
            next;
        }


        ## ------------------------------------------------------------------------------
        ## Classify the PE reads into DEL, INS, WT or other-mismatches. Stores the DEL and INS in
        ## the $data_file
        ## ------------------------------------------------------------------------------
        if ($md1 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $md2 =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $cigar1 =~ /^\d+M\d*D\d+M$/ and $cigar2 =~ /^\d+M\d*D\d+M$/) {
            $stats->{$STAT_OK_DELETION}++;
            my ($deletion_length) = $cigar1 =~ /(\d+)D/;
#             my ($insertion_length2) = $cigar2 =~ /(\d+)D/;
#             my ($insertion_length3) = length(($md1 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#             my ($insertion_length4) = length(($md2 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/)[0]);
#             print join("--", $insertion_length, $insertion_length2, $insertion_length3, $insertion_length4, $cigar1, $cigar2, $md1, $md2), "\n";
#             die if ($insertion_length != $insertion_length2);
#             die if ($insertion_length != $insertion_length3);
#             die if ($insertion_length != $insertion_length4);
            my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
            my ($deletion_seq) = $md1 =~ /^MD:Z:\d+\^([A-Z]+)\d+$/;
            print DATA join("\t", $qname1, "DEL", $deletion_length, $position, $position + $deletion_length, $position + $deletion_length/2, $deletion_seq), "\n";
        } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M\d*I\d+M$/ and $cigar2 =~ /^\d+M\d*I\d+M$/) { 
            $stats->{$STAT_OK_INSERTION}++;
            my ($insertion_length) = $cigar1 =~ /(\d+)I/;
            my $position = ($cigar1 =~ /^(\d+)M/)[0] + $pos1;
            my $insertion_seq = substr($original_seq1, $position, $insertion_length);
            print DATA join("\t", $qname1, "INS", $insertion_length, $position + 1, $position, $position + 1/2, $insertion_seq), "\n";
        } elsif ($md1 =~ /^MD:Z:\d+$/ and $md2 =~ /^MD:Z:\d+$/ and $cigar1 =~ /^\d+M$/ and $cigar2 =~ /^\d+M$/) {
            $stats->{$STAT_OK_WILD_TYPE}++;
            next;
        } else {
#             print join("\t", $cigar1, $md1, $cigar2, $md2), "\n";
            $stats->{$STAT_OTHER_MISMATCHES}++;
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
  Arg[4]        : hashref $stats
  Arg[5]        : num $wt
  Arg[5]        : string $output_pdf_file
  Example       : my $R_file = write_R_script("test", "file.bam.data.txt", $stats, $top_sequences,
                    3123, "test.pdf")
  Description   : Creates an R script to plot all the insertions and deletions in PDF, PNG and SVG
                  format
  Returns       : string with the filename with the resulting R code
  Exceptions    : 

=cut

sub write_R_script {
    my ($label, $data_file, $ref_seq, $top_sequences, $stats, $output_pdf_file) = @_;

    my $num_lines = @$top_sequences;
    my $R_script = $data_file;
    $R_script =~ s/.data.txt$/.R/;

    my $wt = ($stats->{$STAT_OK_WILD_TYPE} or 0);

    open(R, ">$R_script") or die;
    print R "
# =============================================================================
#  Input data
# =============================================================================
#  Data is a matrix of events for the PE reads that pass all the filters.
#  Stats is a matrix with stats re the number of events found in the bam file.
# =============================================================================
data <- read.table('$data_file', header=T, row.names=1)
stats = matrix(data=c(\"",
    join("\",\"", @STAT_ORDER,
        (map {$stats->{$_} or 0} @STAT_ORDER),
        (map {$STAT_COLOR->{$_} or 0} @STAT_ORDER)),
         "\"), ncol=3, byrow=F)
stats = subset(stats, stats[,2]>0);
stats.fail = subset(stats, stats[,1]!='$STAT_OK_WILD_TYPE' & stats[,1]!='$STAT_OK_DELETION' & stats[,1]!='$STAT_OK_INSERTION')
stats.ok = subset(stats, stats[,1]=='$STAT_OK_WILD_TYPE' | stats[,1]=='$STAT_OK_DELETION' | stats[,1]=='$STAT_OK_INSERTION')

force.plots = $force_plots
plot.pdf = $plot_pdf
plot.png = $plot_png
plot.svg = $plot_svg

ref.seq = '$ref_seq'

label = '$label'

# =============================================================================
#  Colors
# =============================================================================
#  Get the colors from the Perl script. Also gets darker versions for
#  insertions and deletions.
# =============================================================================
darken.color <- function(col, amount) {
    mix.col.rgb = colorRamp(c(col, 'black'))(amount)[1,]/255
    mix.col = rgb(mix.col.rgb[1],mix.col.rgb[2],mix.col.rgb[3])
    return(mix.col)
}
col.del = '$COLOR_DELETION'
col.ins = '$COLOR_INSERTION'
dark.col.del = darken.color(col.del, 0.3);
dark.col.ins = darken.color(col.ins, 0.3);

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
num.nhej = num.del + num.ins
num.total = num.wt + num.del + num.ins

perc.wt = paste0(format(100*num.wt/num.total, digits=3),'%')
perc.del = paste0(format(100*num.del/num.total, digits=3),'%')
perc.ins = paste0(format(100*num.ins/num.total, digits=3),'%')
perc.nhej = paste0(format(100*(num.ins+num.del)/num.total, digits=3),'%')


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
    h.del = hist(as.numeric(data.del[,1]), breaks=breaks, plot=F);
    h.ins = hist(as.numeric(data.ins[,1]), breaks=breaks, plot=F);

    ### Cosmetic variables
    ylim.max = max(h.del\$counts, h.ins\$counts)
    xlim = c(min(h.del\$breaks)+1,max(h.del\$breaks))
    xlab.del = paste0('Deletion sizes (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Insertion sizes (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    ### Plot histogram for deletions only (if any)
    main=paste0('Histogram of Deletion sizes (', label, ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.del\$mids, 0, h.del\$mids+1, h.del\$counts, col=col.del)
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    ### Plot histogram for insertions only (if any)
    main=paste0('Histogram of Insertion sizes (', label, ')')
    if (num.ins > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(h.ins\$mids, 0, h.ins\$mids+1, h.ins\$counts, col=col.ins)
    } else if (force.plots) {
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
    mtext(paste0('Histogram of event sizes (', label, ')'), side=3, line = 5, cex=1.2, font=2)
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
plot.location.images <- function(data.del, data.ins, ref.seq='') {

    ### Use same breaks (x-axis) for all 3 plots. Make sure the range spans at least 10 bp.
    if (length(data.del[,4]) + length(data.ins[,4]) > 0 ) {
        ## For upper limit, use data.ins[,4]+1 as insertions mid-points are 0.5 before the insertion
        breaks = as.integer(min(data.del[,4],data.ins[,4], na.rm=T)-1):max(data.del[,4], data.ins[,4]+1, na.rm=T)+1
    } else {
        breaks = 1:10
    }
    while (length(breaks) < 10) {
        if (breaks[1] > 0) {
            breaks = c(breaks[1]-1, breaks)
        }
        breaks = c(breaks, breaks[length(breaks)]+1)
    }

    ### Get the histograms for deletions and insertions, using the set breaks. Store the
    ### result instead of plotting it.
    h.del = hist(as.numeric(data.del[,4]), breaks=breaks, plot=F);
    h.ins = hist(as.numeric(data.ins[,4]), breaks=(breaks-0.5)[1:length(breaks)-1], plot=F);

    ### Cosmetic variables
    ylim.max = 1.04*max(h.del\$counts, h.ins\$counts)
    xlim = c(breaks[1],breaks[length(breaks)])
    xlab.del = paste0('Location of Deletions (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Location of Insertions (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    ### Coordinates are 1-based inclusive. starts array marks the start of the nucleotide. For
    ### instance, bp 100 starts at 100 and finishes at 100+1. We use starts for the deletions
    ### and for the sequence. The insertions happen necessarily in between two nucleotides and
    ### this is why we use the mid-points for these.
    starts = breaks[1]:breaks[length(breaks)-1]

    ### Plot histogram for deletions only (if any)
    main=paste0('Location of the deletions (', label, ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        rect(starts, 0, starts+1, h.del\$counts, col=col.del)
        if (nchar(ref.seq) > 0) {
            cex = 40/max(40, xlim[2]-xlim[1]+1)
            text(starts+0.5, rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]-1), '')[[1]], adj=c(0.5,1.3), cex=cex)
        }
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    ### Plot histogram for insertions only (if any)
    main=paste0('Location of the insertions (', label, ')')
    if (num.ins) {
        plot(h.ins\$counts, xlim=xlim, type='n', xlab=xlab.ins, ylab=ylab, main=main)
        rect(h.ins\$mids+0.5, 0, h.ins\$mids+1.5, h.ins\$counts, col=col.ins)
        if (nchar(ref.seq) > 0) {
            cex = 40/max(40, xlim[2]-xlim[1]+1)
            text(starts+0.5, rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]-1), '')[[1]], adj=c(0.5,1.3), cex=cex)
        }
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No insertions'))
    }

    ### Plot histograms for both Deletion and Insertions (as a mirror plot)
    # Change the margins to make room for the upper axis + label + title
    mar = par('mar')
    par('mar' = c(mar[1]-0.5, mar[2], mar[3]+3, mar[4]))
    # Plot the data
    plot(h.del\$counts, xlim=xlim, type='n', xlab=NA, ylab=ylab,
        main=NA, ylim=c(-ylim.max, ylim.max), yaxt='n');
    if (nchar(ref.seq) > 0) {
        factor = 10/100*40/max(40, xlim[2]-xlim[1]+1)
        usr = par('usr');
        par('usr' = c(usr[1], usr[2], usr[3]*(1+factor), usr[4]))
        ax.ticks = axTicks(2);
        axis(2, col=dark.col.del, col.axis=dark.col.del, at=subset(ax.ticks, ax.ticks>=0))
        rect(starts, 0, starts+1, h.del\$counts, col=col.del)
        cex = 40/max(40, xlim[2]-xlim[1]+1)
        text(starts+0.5, rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]-1), '')[[1]], adj=c(0.5,1.4), cex=cex)
        par('usr' = c(usr[1], usr[2], usr[3], usr[4] *(1+factor)))
        axis(2, col=dark.col.ins, col.axis=dark.col.ins, at=subset(ax.ticks, ax.ticks<=0))
        rect(h.ins\$mids+0.5, 0, h.ins\$mids+1.5, -h.ins\$counts, col=col.ins)
    } else {
        ax.ticks = axTicks(2);
        axis(2, col=dark.col.del, col.axis=dark.col.del, at=subset(ax.ticks, ax.ticks>=0))
        rect(h.del\$mids-1, 0, h.del\$mids, h.del\$counts, col=col.del)
        axis(2, col=dark.col.ins, col.axis=dark.col.ins, at=subset(ax.ticks, ax.ticks<=0))
        rect(h.ins\$breaks, 0, h.ins\$breaks+1, c(-h.ins\$counts,0), col=col.ins)
    }
    # Add the top axis
    axis(3)
    # Write the main title (as an mtext so we can have it a little higher than by default)
    mtext(paste0('Midpoint location (', label, ')'), side=3, line = 5, cex=1.2, font=2)
    # Write the labels for the top and bottom x-axes
    mtext(xlab.del, side=3, line = 2.5)
    mtext(xlab.ins, side=1, line = 2.5)
    # Reset margins to previous value
    par('mar' = mar)


    ### Plot the scatter plot for deletions (if any). Uses the same x-axis as the other plots
    main = paste0('Scatter plot of Deletions sizes vs Midpoint location (', label, ')')
    if (num.del > 0) {
        ylim = c(0,max(data.del[,1],10))
        smoothScatter(data.del[,4], data.del[,1], xlim=xlim, ylim=ylim,
            main=main, xlab=xlab.del, ylab='Deletion size',
            colramp = colorRampPalette(c('white', dark.col.del)));
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }
    
    ### Plot the scatter plot for insertions (if any). Uses the same x-axis as the other plots
    main = paste0('Scatter plot of Insertion sizes vs Midpoint location (', label, ')')
    if (num.ins > 0) {
        ylim = c(0,max(data.ins[,1],10))
        smoothScatter(data.ins[,4], data.ins[,1], xlim=xlim, ylim=ylim,
            main=main, xlab=xlab.ins, ylab='Insertion size',
            colramp = colorRampPalette(c('white', dark.col.ins)));
    } else if (force.plots) {
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
    main=paste0('Frequency of deletion per bp (', label, ')')

    if (nrow(range) == 0) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
        return()
    }

    ### Create a vector to store the number of times each bp has been deleted
    del <- vector(mode='numeric', length=max(data[,3]))
    for (i in 1:dim(range)[1]) { for (p in range[i,1]:range[i,2]) { del[p] <- del[p]+1 } }

    ### Genereate the plot
    # empty plot
    xlim = c(min(range\$from),max(range\$to)+1)
    plot(del, xlim=xlim, type='n',
        main=main,
        xlab='Location of deleted bp', ylab='counts')
    # add light blue rectangles with white background
    rect(min(range\$from):max(range\$to), 0, min(range\$from):max(range\$to)+1,
        del[min(range\$from):max(range\$to)], col=col.del, border='white')
    # add a black line around the profile of deletion frequencies
    lines(c(min(range\$from), min(range\$from):(max(range\$to)+1),min(range\$from)),
        c(0, del[min(range\$from):max(range\$to)], 0, 0), type='s', col='black')
    # Add ref sequence if available
        if (nchar(ref.seq) > 0) {
            cex = 40/max(40, xlim[2]-xlim[1]+1)
            text(xlim[1]:(xlim[2]-1)+0.5, rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]-1), '')[[1]], adj=c(0.5,1.3), cex=cex)
        }
#     if (nchar(ref.seq) > 0) {
#         cex = 40/max(40, xlim[2]-xlim[1]+1)
#         text(xlim[1]:xlim[2], rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]), '')[[1]], adj=c(0.5,1.2),cex=cex)
#     }

}

# =============================================================================
#  FUNCTION plot.pie.chart
# =============================================================================
#  This method plots a pie chart with the number of WT, DEL and INS sequences
# =============================================================================
plot.pie.chart <- function() {
    main = paste0('Summary of events (', label, ')')
    my.stats = stats
    if (sum(as.numeric(my.stats[,2])) > 0) {
        pie(as.numeric(my.stats[,2]), labels=my.stats[,1], main=main, col=my.stats[,3])
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No events'))
    }
    if ((sum(as.numeric(my.stats[,2])) > 0) | force.plots) {
        mtext('5\\'/3\\' mut and indels: Differences with the WT seq before the overlap', side=1, line=0)
        mtext('Reads differ: The reads differ on the overlap', side=1, line=1)
        mtext('Other mismatches: Differences with the WT seq in the overlap or more than 1 indel', side=1, line=2)
    }
    my.stats = stats.ok
    if (sum(as.numeric(my.stats[,2])) > 0) {
        pie(as.numeric(my.stats[,2]), labels=my.stats[,1], main=main, col=my.stats[,3])
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No events'))
    }
    if ((sum(as.numeric(my.stats[,2])) > 0) | force.plots) {
        mtext(paste0('Wild-type: n = ', num.wt, '; ', perc.wt), side=1, line=0, adj=0, at=-1.1)
        mtext(paste0('NHEJ: n = ', num.nhej, '; ', perc.nhej), side=1, line=1, adj=0, at=-1.1)
        mtext(paste0('Deletions: n = ', num.del, '; ', perc.del), side=1, line=0, adj=0, at=0.1)
        mtext(paste0('Insertions: n = ', num.ins, '; ', perc.ins), side=1, line=1, adj=0, at=0.1)
    }
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
    par('mar' = c(mar[1],1,mar[3],0.1))
    ### Initiate plot
    plot(NA,xlim=c(0,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA,
        main=paste0('Most common events (', label, ')'))
    ### Add text
    text(0, 1-1:$num_lines/28, cex=0.3, adj=0, family='mono', labels=c('", join("', '", @$top_sequences), "'))
    par('mar' = mar)
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
";
    if (@$top_sequences) {
        print R "   plot.topseqs(data.del)\n";
    }
    print R "

    ### Check that there are data to be plotted
     if (force.plots | num.del+num.ins > 0) {
    
        plot.size.histograms(data.del, data.ins)
        plot.location.images(data.del, data.ins, ref.seq=ref.seq)
        plot.deletion.frequencies(data.del)

    }
}

# =============================================================================
#  Create and save the plots
# =============================================================================
#  Here the plots are created and saved to a PDF file , then to a set of PNG
#  files and finally to a set of SVG files.
# =============================================================================
if (plot.pdf) {
    pdf('$output_pdf_file')
    plot.figures()
    dev.off()
}

if (plot.png) {
    # Substitute the .pdf extension by %02d.png (to create several files like 01, 02, etc)
    png(sub('.pdf', '.%02d.png', '$output_pdf_file'))
    plot.figures();
    dev.off()
}

if (plot.svg) {
    # Substitute the .pdf extension by %02d.svg (to create several files like 01, 02, etc)
    svg(sub('.pdf', '.%02d.svg', '$output_pdf_file'))
    plot.figures();
    dev.off()
}
";

    close(R);
    
    return($R_script);
}
