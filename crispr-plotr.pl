#! /usr/bin/env perl
use warnings;
use strict;

use Pod::Usage;
use Getopt::Long;
use File::Basename;

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

=item B<--guide-seq> guide_seq.fa>

This is a simple FASTA file (header is optional) with one sequence only (typically a very short one)
corresponding to the guide sequence. It must match perfectly the reference sequence.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

    Example:
    ---------------------------------------
    AGGAGGTCCTA
    ---------------------------------------

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
my $guide_seq_file;
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
my $test = 0;

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
my $STAT_OK_COMPLEX = "COM";


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
    $STAT_OK_COMPLEX,
);

my $COLOR_DELETION = "cyan3";
my $COLOR_INSERTION = "darkorchid3";
my $COLOR_COMPLEX = "#4d80cd";
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
    "$STAT_OK_COMPLEX" => $COLOR_COMPLEX,
};


GetOptions(
    "help"  => \$help,
    "debug"  => \$debug,
    "test"  => \$test,
    "label=s" => \$label,
    "input_file|input-file=s" => \$input_bam_file,
    "output_file|output-file=s" => \$output_pdf_file,
    "ref_seq|ref-seq|wt_file|wt-file|wt_seq|wt-seq=s" => \$ref_seq_file,
    "guide_seq|guide-seq=s" => \$guide_seq_file,
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

if ($test) {
    use Test::More;
    test_merge_reads();
    done_testing();
    exit(0);
}

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
my $guide_seq = "";
if ($ref_seq_file) {
    $ref_seq = read_fasta_seq($ref_seq_file);
    
    if ($guide_seq_file) {
        $guide_seq = read_fasta_seq($guide_seq_file);
    }

    my $wt = ($stats->{$STAT_OK_WILD_TYPE} or 0);
    $top_sequences = get_top_sequences($ref_seq, $guide_seq, $data_file, $wt);
}

my $R_script = write_R_script($label, $data_file, $ref_seq, $guide_seq, $top_sequences, $stats, $output_pdf_file);

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

    open(FASTA, $fasta_file) or die "Cannot open FASTA file <$fasta_file>\n";
    while (<FASTA>) {
        chomp;
        next if (/^>/);
        $seq .= uc($_);
    }

    return $seq;
}

=head2 get_top_sequences

  Arg[1]        : string $ref_seq
  Arg[2]        : string $guide_seq
  Arg[3]        : string $data_filename
  Arg[4]        : integer $num_of_wild_type_seqs
  Example       : my $top_sequences = get_top_sequences($ref_seq, $guide_seq, $data_file, 1213);
  Description   : Reads from the $data_file the most common deletions and insertions (up to 10) and
                  aligns them to the wild-type sequence. Highlights the guide seq if provided.
  Returns       : arrayref of strings
  Exceptions    : prints WARNING if $guide_seq does not match $ref_seq

=cut

sub get_top_sequences {
    my ($ref_seq, $guide_seq, $data_file, $wt) = @_;
    my $top_sequences = [];


    ## ------------------------------------------------------------------------------
    ## Highlight the guide sequence (if provided) as uppercase vs lowercase
    ## ------------------------------------------------------------------------------
    my $guide_start;
    if ($guide_seq) {
        $guide_start = index($ref_seq, $guide_seq);
        if ($guide_start >= 0) {
            substr($ref_seq, 0, $guide_start) = lc(substr($ref_seq, 0, $guide_start));
            substr($ref_seq, $guide_start+length($guide_seq)) = lc(substr($ref_seq, $guide_start+length($guide_seq)));
            print $ref_seq,"\n";
        } else {
            my $revcom_guide_seq = revcom($guide_seq);
            $guide_start = index($ref_seq, $revcom_guide_seq);
            if ($guide_start >= 0) {
                substr($ref_seq, 0, $guide_start) = lc(substr($ref_seq, 0, $guide_start));
                substr($ref_seq, $guide_start+length($guide_seq)) = lc(substr($ref_seq, $guide_start+length($guide_seq)));
                print $ref_seq,"\n";
            } else {
                print STDERR "======================================================================\n";
                print STDERR "WARNING: guide sequence ('$guide_seq') not found in amplicon sequence\n";
                print STDERR "WARNING: guide sequence ('$revcom_guide_seq') not found in amplicon sequence\n";
                print STDERR "======================================================================\n";
            }
        }
    }

    ## ------------------------------------------------------------------------------
    ## Reads and extract stats on most common deletions and insertions:
    ## ------------------------------------------------------------------------------
    my @del_lines = qx"more $data_file | awk '\$2 == \"DEL\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn";
    my @ins_lines = qx"more $data_file | awk '\$2 == \"INS\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn";
    my @com_lines = qx"more $data_file | awk '\$2 == \"COM\" { print \$3, \$4, \$7}' | sort | uniq -c | sort -rn";


    ## ------------------------------------------------------------------------------
    ## Calculates the region of the REF sequence to display
    ## ------------------------------------------------------------------------------
    my $min_from;
    my $max_to;
    my $longest_seq = 0;
    # Include the guide sequence if available
    if ($guide_start) {
        $min_from = $guide_start + 1;
        $max_to = $guide_start + length($guide_seq);
    }
    # Expand as necessary to cover all deletion and complex cases
    foreach my $this_line (@del_lines, @com_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        $min_from = $from if (!$min_from or $from < $min_from);
        $max_to = $from+$del_length if (!$max_to or $from+$del_length > $max_to);
        $longest_seq = length($seq) if ($longest_seq < length($seq));
    }
    # Same for insertions (except that insertions are of length 0 in ref coordinates)
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        $min_from = $from - 1 if (!$min_from or $from - 1 < $min_from);
        $max_to = $from if (!$max_to or $from > $max_to);
        $longest_seq = length($seq) + 2 if ($longest_seq < length($seq) + 2);
    }
    # Just in case there are no guide sequence nor deletions nor insertions
    if (!defined($max_to)) {
        $max_to = int(length($ref_seq) / 2) + 10;
        $min_from = $max_to - 10;
    }
    # Expand the sequence a little further to give more context
    $min_from -= 20;
    $min_from = 1 if ($min_from < 1);
    $max_to += 19;


    ## ------------------------------------------------------------------------------
    ## Sets the output format (using whitespaces to be nicely printed in R afterewards)
    ## ------------------------------------------------------------------------------
    my $format = "\%-".($max_to - $min_from + 1)."s \%7s  \%-4s \%3s \%6s \%-${longest_seq}s";


    ## ------------------------------------------------------------------------------
    ## Header and WT sequence
    ## ------------------------------------------------------------------------------
    my $header = sprintf($format, "Sequence" , "Num", "TYPE", "L", "POS", "Diff");
    my $wt_sequence = sprintf($format, substr($ref_seq, $min_from - 1, $max_to - $min_from), $wt, "WT", 0, "NA", "");
    $top_sequences = [$header, "", $wt_sequence, ""];


    ## ------------------------------------------------------------------------------
    ## Most common deletions
    ## ------------------------------------------------------------------------------
    my $del_sequences_hash = {};
    foreach my $this_line (@del_lines) {
        my ($num, $del_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, ($min_from - 1), ($from - $min_from)) .
                          '-' x $del_length .
                          substr($ref_seq, ($from + $del_length - 1), $max_to - ($from + $del_length - 1));
        my $resulting_seq = substr($ref_seq, ($min_from - 1), ($from - $min_from)) .
                            substr($ref_seq, ($from + $del_length - 1), $max_to - ($from + $del_length - 1));
        if ($del_sequences_hash->{$resulting_seq}) {
            $del_sequences_hash->{$resulting_seq}->{num} += $num;
        } else {
            $del_sequences_hash->{$resulting_seq}->{aligned_seq} = $aligned_seq;
            $del_sequences_hash->{$resulting_seq}->{num} = $num;
            $del_sequences_hash->{$resulting_seq}->{del_length} = $del_length;
            $del_sequences_hash->{$resulting_seq}->{from} = $from;
            $del_sequences_hash->{$resulting_seq}->{seq} = $seq;
        }
    }
    my @del_sequences = map {
            sprintf($format, $del_sequences_hash->{$_}->{aligned_seq}, $del_sequences_hash->{$_}->{num},
                    "DEL", $del_sequences_hash->{$_}->{del_length}, $del_sequences_hash->{$_}->{from},
                    ">".$del_sequences_hash->{$_}->{seq}."<")
            } (sort {$del_sequences_hash->{$b}->{num} <=> $del_sequences_hash->{$a}->{num} ||
                     $del_sequences_hash->{$a}->{del_length} <=> $del_sequences_hash->{$b}->{del_length}
                    } keys $del_sequences_hash);
    my $deletions_file = $data_file;
    $deletions_file =~ s/\.txt$/.del.txt/;
    open(DEL, ">$deletions_file") or die;
    print DEL join("\n", $header, $wt_sequence, "", @del_sequences, "");
    close(DEL);
    push(@$top_sequences, splice(@del_sequences, 0, 18));


    ## ------------------------------------------------------------------------------
    ## Separation line
    ## ------------------------------------------------------------------------------
    push(@$top_sequences, "");


    ## ------------------------------------------------------------------------------
    ## Most common insertions
    ## ------------------------------------------------------------------------------
    my $ins_sequences_hash = {};
    foreach my $this_line (@ins_lines) {
        my ($num, $ins_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from - 1, $max_to - $min_from + 1);
        substr($aligned_seq, $from - $min_from - 1, 2, "><");
        my $resulting_seq = substr($ref_seq, $min_from - 1, $max_to - $min_from + 1);
        substr($resulting_seq, $from - $min_from, 0, $seq);
        if ($ins_sequences_hash->{$resulting_seq}) {
            $ins_sequences_hash->{$resulting_seq}->{num} += $num;
        } else {
            $ins_sequences_hash->{$resulting_seq}->{aligned_seq} = $aligned_seq;
            $ins_sequences_hash->{$resulting_seq}->{num} = $num;
            $ins_sequences_hash->{$resulting_seq}->{ins_length} = $ins_length;
            $ins_sequences_hash->{$resulting_seq}->{from} = $from;
            $ins_sequences_hash->{$resulting_seq}->{seq} = $seq;
        }
    }
    my @ins_sequences = map {
            sprintf($format, $ins_sequences_hash->{$_}->{aligned_seq}, $ins_sequences_hash->{$_}->{num},
                    "INS", $ins_sequences_hash->{$_}->{ins_length}, $ins_sequences_hash->{$_}->{from},
                    ">".$ins_sequences_hash->{$_}->{seq}."<")
            } (sort {$ins_sequences_hash->{$b}->{num} <=> $ins_sequences_hash->{$a}->{num} ||
                     $ins_sequences_hash->{$a}->{ins_length} <=> $ins_sequences_hash->{$b}->{ins_length}
                    } keys $ins_sequences_hash);
    my $insertions_file = $data_file;
    $insertions_file =~ s/\.txt$/.ins.txt/;
    open(INS, ">$insertions_file") or die;
    print INS join("\n", $header, $wt_sequence, "", @ins_sequences, "");
    close(INS);
    push(@$top_sequences, splice(@ins_sequences, 0, 18));

    ## ------------------------------------------------------------------------------
    ## Separation line
    ## ------------------------------------------------------------------------------
    push(@$top_sequences, "");


    ## ------------------------------------------------------------------------------
    ## Most common complex cases
    ## ------------------------------------------------------------------------------
    my $com_sequences_hash = {};
    foreach my $this_line (@com_lines) {
        my ($num, $com_length, $from, $seq) = $this_line =~ /(\d+)\s(\-?\d+)\s(\d+)\s(\w+)/;
        my $aligned_seq = substr($ref_seq, $min_from - 1, $max_to - $min_from + 1);
        substr($aligned_seq, $from - $min_from - 1, $com_length + 2, ">" . ("-" x $com_length) . "<");
        my $resulting_seq = substr($ref_seq, $min_from - 1, $max_to - $min_from + 1);
        substr($resulting_seq, $from - $min_from, $com_length, $seq);
        if ($com_sequences_hash->{$resulting_seq}) {
            $com_sequences_hash->{$resulting_seq}->{num} += $num;
        } else {
            $com_sequences_hash->{$resulting_seq}->{aligned_seq} = $aligned_seq;
            $com_sequences_hash->{$resulting_seq}->{num} = $num;
            $com_sequences_hash->{$resulting_seq}->{com_length} = $com_length;
            $com_sequences_hash->{$resulting_seq}->{from} = $from;
            $com_sequences_hash->{$resulting_seq}->{seq} = $seq;
        }
    }
    my @com_sequences = map {
            sprintf($format, $com_sequences_hash->{$_}->{aligned_seq}, $com_sequences_hash->{$_}->{num},
                    "COM", $com_sequences_hash->{$_}->{com_length}, $com_sequences_hash->{$_}->{from},
                    ">".$com_sequences_hash->{$_}->{seq}."<")
            } (sort {$com_sequences_hash->{$b}->{num} <=> $com_sequences_hash->{$a}->{num} ||
                     $com_sequences_hash->{$a}->{com_length} <=> $com_sequences_hash->{$b}->{com_length}
                    } keys $com_sequences_hash);

    my $complex_file = $data_file;
    $complex_file =~ s/\.txt$/.com.txt/;
    open(COM, ">$complex_file") or die;
    print COM join("\n", $header, $wt_sequence, "", @com_sequences, "");
    close(COM);
    push(@$top_sequences, splice(@com_sequences, 0, 18));

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
        
        my $read1 = {
            'seq' => $seq1,
            'qual' => $qual1,
            'start' => $pos1,
            'cigar' => $cigar1,
            'md' => $md1};

        my $read2 = {
            'seq' => $seq2,
            'qual' => $qual2,
            'start' => $pos2,
            'cigar' => $cigar2,
            'md' => $md2};

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
        ## Merge the reads and check for exceptions
        ## ------------------------------------------------------------------------------
        my $merged_read = merge_reads($read1, $read2);

        if (exists($merged_read->{'error'})) {
            $stats->{$merged_read->{'error'}}++;
            next;
        }


        if (!$allow_indels and $merged_read->{'5prime'}{'cigar'} !~ /^\d+M$/) {
            $stats->{$STAT_INDEL_5_PRIME}++;
            next;
        }

        if (!$allow_muts and $merged_read->{'5prime'}{'md'} !~ /^MD:Z:\d+$/) {
            $stats->{$STAT_MUTATION_5_PRIME}++;
            next;
        }

        if (!$allow_indels and $merged_read->{'3prime'}{'cigar'} !~ /^\d+M$/) {
            $stats->{$STAT_INDEL_3_PRIME}++;
            next;
        }

        if (!$allow_muts and $merged_read->{'3prime'}{'md'} !~ /^MD:Z:\d+$/) {
            $stats->{$STAT_MUTATION_3_PRIME}++;
            next;
        }

        my $cigar_overlap = $merged_read->{'overlap'}{'cigar'};
        my $md_overlap = $merged_read->{'overlap'}{'md'};
        if ($md_overlap =~ /^MD:Z:\d+\^[A-Z]+\d+$/ and $cigar_overlap =~ /^(\d+M)?\d+D(\d+M)?$/) {
            $stats->{$STAT_OK_DELETION}++;
            my ($deletion_start) = ($cigar_overlap =~ /^(\d+)M/);
            my ($deletion_length) = $cigar_overlap =~ /(\d+)D/;
            my ($deletion_seq) = $md_overlap =~ /^MD:Z:\d+\^([A-Z]+)\d+$/;
            my $from = $merged_read->{'overlap'}{'start'} + ($deletion_start or 0);
            my $to = $from + $deletion_length - 1;
            my $midpoint = $from + ($deletion_length - 1)/2;
            print DATA join("\t", $qname1, "DEL", $deletion_length, $from, $to, $midpoint, $deletion_seq), "\n";
        } elsif ($md_overlap =~ /^MD:Z:\d+$/ and $cigar_overlap =~ /^\d+M\d*I\d+M$/) {
            $stats->{$STAT_OK_INSERTION}++;
            my ($insertion_start) = $cigar_overlap =~ /^(\d+)M/;
            my ($insertion_length) = $cigar_overlap =~ /(\d+)I/;
            my $from = $merged_read->{'overlap'}{'start'} + ($insertion_start or 0);
            my $to = $from - 1;
            my $midpoint = $from - 1/2;
            my $insertion_seq = substr($merged_read->{'overlap'}{'seq'}, ($insertion_start or 0), $insertion_length);
            print DATA join("\t", $qname1, "INS", $insertion_length, $from, $to, $midpoint, $insertion_seq), "\n";
        } elsif ($md_overlap =~ /^MD:Z:\d+$/ and $cigar_overlap =~ /^\d+M$/) {
            $stats->{$STAT_OK_WILD_TYPE}++;
        } elsif ($cigar_overlap !~ /^\d+M$/) {
            $stats->{$STAT_OK_COMPLEX}++;
            my ($event_5prime_md) = ($md_overlap =~ /^MD:Z:(\d+)/);
            $event_5prime_md ||= 0;
            my ($event_5prime_cigar) = ($cigar_overlap =~ /^(\d+)M/);
            $event_5prime_cigar ||= 0;
            # $event_5prime_length is the length of the sequence in the overlap that is identical to the WT seq
            my $event_5prime_length = $event_5prime_md>$event_5prime_cigar?$event_5prime_cigar:$event_5prime_md;

            my ($event_3prime_md) = ($md_overlap =~ /(\d+)$/);
            $event_3prime_md ||= 0;
            my ($event_3prime_cigar) = ($cigar_overlap =~ /(\d+)M$/);
            $event_3prime_cigar ||= 0;
            # $event_3prime_length is the length of the sequence in the overlap that is identical to the WT seq
            my $event_3prime_length = $event_3prime_md>$event_3prime_cigar?$event_3prime_cigar:$event_3prime_md;

            my $from = $merged_read->{'overlap'}{'start'} + $event_5prime_length;
            my $to = $merged_read->{'overlap'}{'end'} - $event_3prime_length;
            my $midpoint = ($from + $to)/2;
            my $complex_seq = $merged_read->{'overlap'}{'seq'};
            substr($complex_seq, 0, $event_5prime_length, "");
            substr($complex_seq, -$event_3prime_length, $event_3prime_length, "");
            print DATA join("\t", $qname1, "COM", ($to - $from + 1), $from, $to, $midpoint, $complex_seq), "\n";
#             print "$cigar_overlap $md_overlap  $event_5prime_length-$event_3prime_length $complex_seq\n";
#             <STDIN>;
        }
    }
    close(SAM);
    close(DATA);

    return ($data_file, $stats);
}

sub revcom {
    my ($seq) = @_;
    $seq = reverse($seq);
    $seq =~ tr/ACTG/TGAC/;
    return $seq;
}

sub get_cigar_expanded_string {
    my ($cigar_compact_string) = @_;
    my $cigar_expanded_string = "";

    my @cigars = ($cigar_compact_string =~ /(\d*)(\w)/g);
    for (my $i=0; $i < @cigars; $i+=2) {
        $cigar_expanded_string .= $cigars[$i+1] x $cigars[$i];
    }

    return $cigar_expanded_string;
}

sub get_cigar_compact_string {
    my ($cigar_expanded_string) = @_;
    my $cigar_compact_string = "";

    while ($cigar_expanded_string =~ /^((.)\g2*)/) {
        $cigar_compact_string .= length($1).substr($1,0,1);
        $cigar_expanded_string =~ s/^$1//;
    }

    return $cigar_compact_string;
}

sub get_md_expanded_string {
    my ($md_compact_string) = @_;
    my $md_expanded_string = "";

    $md_compact_string =~ s/^MD:Z://;
    my @mds = grep {$_} split(/(\d+)/, $md_compact_string);
    foreach my $md_bit (@mds) {
        if ($md_bit =~ /^\d+$/) {
            # Matches: write dots
            $md_expanded_string .= "." x $md_bit;
        } elsif ($md_bit =~ /^[A-Z]+$/) {
            # Mismatches: write lowercase
            $md_expanded_string .= lc($md_bit);
        } elsif ($md_bit =~ /^\^[A-Z]+$/) {
            # Deletions: write uppercase
            $md_expanded_string .= substr($md_bit, 1);
        } else {
            die "Cannot understand MD:Z flag.\n";
        }
    }
    $md_expanded_string =~ s/n/./g;

    return $md_expanded_string;
}

sub get_md_compact_string {
    my ($md_expanded_string) = @_;
    my $md_compact_string = "MD:Z:";

    my $last_md_bit_mode = "";
    while ($md_expanded_string =~ /^(\^?(.)\g2*)/) {
        my $md_bit = $1;
        if ($md_bit =~ /^\.+$/) {
            $md_compact_string .= length($md_bit);
            $last_md_bit_mode = "match";
        } elsif ($md_bit =~ /^[a-z]+$/) {
            if ($last_md_bit_mode eq "deletion") {
                $md_compact_string .= "0";
            }
            $md_compact_string .= join("0", split("", uc($md_bit)));
            $last_md_bit_mode = "mismatch";
        } elsif ($md_bit =~ /^[A-Z]+$/) {
            $md_compact_string .= ($last_md_bit_mode ne "deletion"?"^":"").$md_bit;
            $last_md_bit_mode = "deletion";
        } else {
            die "Cannot understand MD:Z flag: $md_expanded_string.\n";
        }
        $md_expanded_string =~ s/^$md_bit//;
    }

    return $md_compact_string;
}

sub merge_reads {
    my ($read1, $read2) = @_;

    ## ------------------------------------------------------------------------
    ## Extract values from structure
    ## ------------------------------------------------------------------------
    my $start1 = $read1->{start};
    my $cigar1 = $read1->{cigar};
    my $md1 = $read1->{md};
    my $seq1 = $read1->{seq};
    my $qual1 = $read1->{qual};

    my $start2 = $read2->{start};
    my $cigar2 = $read2->{cigar};
    my $md2 = $read2->{md};
    my $seq2 = $read2->{seq};
    my $qual2 = $read2->{qual};


    ## ------------------------------------------------------------------------
    ## Swap read1 and read2 if necessary
    ## ------------------------------------------------------------------------
    if ($start2 < $start1) {
        ($start1, $start2) = ($start2, $start1);
        ($cigar1, $cigar2) = ($cigar2, $cigar1);
        ($md1, $md2) = ($md2, $md1);
        ($seq1, $seq2) = ($seq2, $seq1);
        ($qual1, $qual2) = ($qual2, $qual1);
    }


    ## ------------------------------------------------------------------------
    ## Get expanded cigar and md strings
    ## ------------------------------------------------------------------------
    my $cigar_str1 = get_cigar_expanded_string($cigar1);
    my $cigar_str2 = get_cigar_expanded_string($cigar2);

    my $md_str1 = get_md_expanded_string($md1);
    my $md_str2 = get_md_expanded_string($md2);


    ## ------------------------------------------------------------------------
    ## Get some coordinates for unique and overlapping sequence
    ## ------------------------------------------------------------------------
    my $end1 = $start1 + length($seq1) - ($cigar_str1 =~ tr/I/I/) + ($cigar_str1 =~ tr/D/D/) - 1;
    my $end2 = $start2 + length($seq2) - ($cigar_str2 =~ tr/I/I/) + ($cigar_str2 =~ tr/D/D/) - 1;
    my $ref_length_5prime_unique = $start2 - $start1;
    my $ref_length_3prime_unique = $end2 - $end1;
    my $ref_start_overlap = $start2;
    my $ref_end_overlap = $end1;
    my $ref_length_overlap = $end1 - $start2 + 1;

    if ($debug) {
        print "\n\n###########################################################################\n";
        print join("\n", "$start1-$end1 in ref coordinates", $seq1, $qual1, $cigar_str1." ($cigar1)", $md_str1." ($md1)"), "\n";
        print "###########################################################################\n";
        print join("\n", "$start2-$end2 in ref coordinates", $seq2, $qual2, $cigar_str2." ($cigar2)", $md_str2." ($md2)"), "\n";
        print "###########################################################################\n\n";
    }

    if ($end1 - $start2 + 1 < $min_overlap) {
        # Reads don't overlap.
        return {'error' => $STAT_NO_OVERLAP};
    }

    if ($end1 > $end2) {
        # Read2 fully included in read1: skip this odd pair.
        return {'error' => $STAT_NO_OVERLAP};
    }


    ## ------------------------------------------------------------------------
    ## Get info for the 5' unique sequence
    ## ------------------------------------------------------------------------
    my ($cigar_5prime_unique) = ($cigar_str1 =~ /^((?:I*[MD]){$ref_length_5prime_unique}I*)/);
    my $seq_length_5prime_unique = ($cigar_5prime_unique =~ tr/MI/MI/);

    ## If overlap starts in an insertion:
    my $length_insertion_end_5prime_unique = length(($cigar_5prime_unique =~ /(I*)$/)[0]);
    my $length_insertion_start_read2 = length(($cigar_str2 =~ /^(I*)/)[0]);
    if ($length_insertion_end_5prime_unique > 0 and $length_insertion_start_read2 > 0) {
        if ($debug) {
            print " -- Trimming end of 5' unique required: ${length_insertion_end_5prime_unique}I vs ${length_insertion_start_read2}I\n";
        }
        if ($length_insertion_end_5prime_unique >= $length_insertion_start_read2) {
            # Insertion in 5' unique is the same or longer: just trim it  to remove the overlapping sequence
            $cigar_5prime_unique =~ s/I{$length_insertion_start_read2}$//;
            $seq_length_5prime_unique -= $length_insertion_start_read2;
        } else {
            return {'error' => $STAT_READ_MISMATCH};
        }
    }
    my ($seq_5prime_unique) = substr($seq1, 0, $seq_length_5prime_unique);
    my ($qual_5prime_unique) = substr($qual1, 0, $seq_length_5prime_unique);
    my $md_length_5prime_unique = ($cigar_5prime_unique =~ tr/MD/MD/);
    my ($md_5prime_unique) = substr($md_str1, 0, $md_length_5prime_unique);

    if ($debug) {
        print "5' sequence:\n$seq_5prime_unique\n$qual_5prime_unique\n$cigar_5prime_unique\n$md_5prime_unique\n\n";
    }


    ## ------------------------------------------------------------------------
    ## Get info for the 3' unique sequence
    ## ------------------------------------------------------------------------
    my ($cigar_3prime_unique) = ($cigar_str2 =~ /(I*(?:I*[MD]){$ref_length_3prime_unique})$/);
    my $seq_length_3prime_unique = ($cigar_3prime_unique =~ tr/MI/MI/);

    ## If overlap ends in an insertion:
    my $length_insertion_start_3prime_unique = length(($cigar_3prime_unique =~ /^(I*)/)[0]);
    my $length_insertion_end_read1 = length(($cigar_str1 =~ /(I*)$/)[0]);
    if ($length_insertion_start_3prime_unique > 0 and $length_insertion_end_read1 > 0) {
        if ($debug) {
            print " -- Trimming start of 3' unique required: ${length_insertion_start_3prime_unique}I vs ${length_insertion_end_read1}I\n";
        }
        if ($length_insertion_start_3prime_unique >= $length_insertion_end_read1) {
            # Insertion in 3' unique is the same or longer: just trim it to remove the overlapping sequence
            $cigar_3prime_unique =~ s/^I{$length_insertion_end_read1}//;
            $seq_length_3prime_unique -= $length_insertion_end_read1;
        } else {
            return {'error' => $STAT_READ_MISMATCH};
        }
    }
    my ($seq_3prime_unique) = substr($seq2, -$seq_length_3prime_unique);
    my ($qual_3prime_unique) = substr($qual2, -$seq_length_3prime_unique);
    my $md_length_3prime_unique = ($cigar_3prime_unique =~ tr/MD/MD/);
    my ($md_3prime_unique) = substr($md_str2, -$md_length_3prime_unique);

    if ($debug) {
        print "3' sequence:\n$seq_3prime_unique\n$qual_3prime_unique\n$cigar_3prime_unique\n$md_3prime_unique\n\n";
    }


    ## ------------------------------------------------------------------------
    ## Get info for the overlapping sequence
    ## ------------------------------------------------------------------------
    my $seq_length_overlapping_sequence1 = length($seq1) - $seq_length_5prime_unique;
    my $seq_length_overlapping_sequence2 = length($seq2) - $seq_length_3prime_unique;
    my $seq_overlap1 = substr($seq1, $seq_length_5prime_unique);
    my $seq_overlap2 = substr($seq2, 0, -$seq_length_3prime_unique);
    my $qual_overlap1 = substr($qual1, $seq_length_5prime_unique);
    my $qual_overlap2 = substr($qual2, 0, -$seq_length_3prime_unique);
    my $cigar_overlap1 = $cigar_str1;
    $cigar_overlap1 =~ s/^$cigar_5prime_unique//;
    my $cigar_overlap2 = $cigar_str2;
    $cigar_overlap2 =~ s/$cigar_3prime_unique$//;
    my $md_overlap1 = $md_str1;
    $md_overlap1 =~ s/^$md_5prime_unique//;
    my $md_overlap2 = $md_str2;
    $md_overlap2 =~ s/$md_3prime_unique$//;
    if ($debug) {
        print "overlapping sequence:\n$seq_overlap1\n$seq_overlap2\n",
                "$qual_overlap1\n$qual_overlap2\n",
                "$cigar_overlap1\n$cigar_overlap2\n",
                "$md_overlap1\n$md_overlap2\n",
                ;
    }


#     if ($cigar_overlap1 ne $cigar_overlap2) {
#         print "Cigar strings don't match\n";
#     }
#     if ($md_overlap1 ne $md_overlap2) {
#         print "MD strings don't match\n";
#     }
#     if ($seq_length_overlapping_sequence1 != $seq_length_overlapping_sequence2) {
#         print "Error parsing sequences. Overlapping sequence lengths do not match\n";
#         <STDIN>;
#         return {'error' => "Overlapping sequences do not match"};
#     }
    if ($debug and 
        ($cigar_overlap1 ne $cigar_overlap2)
        and
        ($seq_overlap1 eq $seq_overlap2)
        ) {
        print "Cigar strings do not match even if sequences do";
        <STDIN>;
    }
    if ($seq_overlap1 ne $seq_overlap2) {
#         print "Sequences don't match\n";
        return {'error' => $STAT_READ_MISMATCH};
    }
    my $md_pos = 0;


    ## ------------------------------------------------------------------------
    ## Get the merged information
    ## ------------------------------------------------------------------------
    my $merged_sequence = undef;
    my $merged_cigar = undef;
    my $merged_md = undef;
    my $merged_qual = undef;
    my $qual_overlap_max = undef;
    if ($seq_overlap1 eq $seq_overlap2) {
        $merged_cigar = get_cigar_compact_string($cigar_5prime_unique.$cigar_overlap1.$cigar_3prime_unique);
        $merged_md = get_md_compact_string($md_5prime_unique.$md_overlap1.$md_3prime_unique);
        $merged_sequence = $seq_5prime_unique.$seq_overlap1.$seq_3prime_unique;
        $qual_overlap_max = "";
        for (my $a = 0; $a < $seq_length_overlapping_sequence1; $a++) {
            my $bp1 = substr($seq_overlap1, $a, 1);
            my $bp2 = substr($seq_overlap2, $a, 1);
            my $q1 = substr($qual_overlap1, $a, 1);
            my $q2 = substr($qual_overlap2, $a, 1);
            $qual_overlap_max .= ord($q1)>ord($q2)?$q1:$q2;
    #         print STDERR join("\t", $bp1, $bp2, ord($q1)-33, ord($q2)-33), "\n";
        }
        $merged_qual = $qual_5prime_unique.$qual_overlap_max.$qual_3prime_unique;
    }

    return {
            'merged' => {
                'start' => $start1,
                'end'   => $end2,
                'seq'   => $merged_sequence,
                'qual'  => $merged_qual,
                'cigar' => $merged_cigar,
                'md'    => $merged_md,
            },
            '5prime' => {
                'start' => $start1,
                'end'   => $ref_start_overlap - 1,
                'seq'   => $seq_5prime_unique,
                'qual'  => $qual_5prime_unique,
                'cigar' => get_cigar_compact_string($cigar_5prime_unique),
                'md'    => get_md_compact_string($md_5prime_unique),
            },
            'overlap' => {
                'start' => $ref_start_overlap,
                'end'   => $ref_end_overlap,
                'seq'   => $seq_overlap1,
                'qual'  => $qual_overlap_max,
                'cigar' => get_cigar_compact_string($cigar_overlap1),
                'md'    => get_md_compact_string($md_overlap1),
            },
            '3prime' => {
                'start' => $ref_end_overlap + 1,
                'end'   => $end2,
                'seq'   => $seq_3prime_unique,
                'qual'  => $qual_3prime_unique,
                'cigar' => get_cigar_compact_string($cigar_3prime_unique),
                'md'    => get_md_compact_string($md_3prime_unique),
            },
        
        };
}


=head2 write_R_script

  Arg[1]        : string $label
  Arg[2]        : string $data_filename
  Arg[3]        : string $ref_seq
  Arg[4]        : string $guide_seq
  Arg[5]        : arrayref $top_sequences
  Arg[6]        : hashref $stats
  Arg[7]        : string $output_pdf_file
  Example       : my $R_file = write_R_script("test", "file.bam.data.txt", "ACTGTGCATGGAATTGGGAACC",
                    "GAATTGG", $top_sequences, $stats, "test.pdf")
  Description   : Creates an R script to plot all the insertions and deletions in PDF, PNG and SVG
                  format
  Returns       : string with the filename with the resulting R code
  Exceptions    : 

=cut

sub write_R_script {
    my ($label, $data_file, $ref_seq, $guide_seq, $top_sequences, $stats, $output_pdf_file) = @_;

    my $num_lines = @$top_sequences;
    my $R_script = $data_file;
    $R_script =~ s/.data.txt$/.R/;

    my $wt = ($stats->{$STAT_OK_WILD_TYPE} or 0);

    my $revcom_guide_seq = revcom($guide_seq);

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
stats.fail = subset(stats, stats[,1]!='$STAT_OK_WILD_TYPE' & stats[,1]!='$STAT_OK_DELETION' & stats[,1]!='$STAT_OK_INSERTION' & stats[,1]!='$STAT_OK_COMPLEX')
stats.ok = subset(stats, stats[,1]=='$STAT_OK_WILD_TYPE' | stats[,1]=='$STAT_OK_DELETION' | stats[,1]=='$STAT_OK_INSERTION' | stats[,1]=='$STAT_OK_COMPLEX')

force.plots = $force_plots
plot.pdf = $plot_pdf
plot.png = $plot_png
plot.svg = $plot_svg

# =============================================================================
#  Sequences
# =============================================================================
ref.seq = '$ref_seq'
guide.seq = '$guide_seq'
guide.seq.revcom = '$revcom_guide_seq'
guide.start = regexpr(guide.seq, ref.seq, fixed=T)[1];
if (guide.start == -1) {
    guide.start = regexpr(guide.seq.revcom, ref.seq, fixed=T)[1];
}
guide.length = nchar(guide.seq)

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
col.com = '$COLOR_COMPLEX'
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
data.com <- subset(data, event=='COM', select=2:ncol(data))


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
num.com = dim(data.com)[1]
num.nhej = num.del + num.ins + num.com
num.total = num.wt + num.del + num.ins + num.com

perc.wt = paste0(format(100*num.wt/num.total, digits=3),'%')
perc.del = paste0(format(100*num.del/num.total, digits=3),'%')
perc.ins = paste0(format(100*num.ins/num.total, digits=3),'%')
perc.com = paste0(format(100*num.com/num.total, digits=3),'%')
perc.nhej = paste0(format(100*(num.nhej)/num.total, digits=3),'%')


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
        plot(h.ins\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
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
        breaks = as.integer(min(data.del[,4],data.ins[,4], na.rm=T)):(max(data.del[,4], data.ins[,4]+1, na.rm=T)+1)
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
    h.del = hist(as.numeric(data.del[,4]), breaks=breaks, right=F, plot=F);
    ### All mid-points for insertions are n+0.5 by definition:
    breaks.ins = (breaks+0.5)[1:length(breaks)-1]
    h.ins = hist(as.numeric(data.ins[,4]), breaks=breaks.ins, right=F, plot=F);

    ### Cosmetic variables
    ylim.max = 1.04*max(h.del\$counts, h.ins\$counts)
    xlim = c(breaks[1],breaks[length(breaks)]-1)
    xlab.del = paste0('Location of Deletions (n = ', num.del, '; ', perc.del, ')')
    xlab.ins = paste0('Location of Insertions (n = ', num.ins, '; ', perc.ins, ')')
    ylab='counts'

    ### Coordinates are 1-based inclusive. starts array marks the start of the nucleotide. For
    ### instance, bp 100 starts at 100 and finishes at 100+1. We use starts for the deletions
    ### and for the sequence. The insertions happen necessarily in between two nucleotides and
    ### this is why we use the mid-points for these.
    starts = breaks[1]:breaks[length(breaks)-1]
    starts.ins = breaks.ins[1]:breaks.ins[length(breaks.ins)-1]

    ### Plot histogram for deletions only (if any)
    main=paste0('Location of the deletions (', label, ')')
    if (num.del > 0) {
        plot(h.del\$counts, xlim=xlim, type='n', xlab=xlab.del, ylab=ylab, main=main)
        if (nchar(ref.seq) > 0) {
            cex = 40/max(40, xlim[2]-xlim[1]+1)
            if (guide.start > 0 && guide.length > 0) {
                rect(starts[1] + guide.start - xlim[1] - 0.5, par('usr')[3]*0.975, starts[1] + guide.start - xlim[1] + guide.length - 0.5, 0, col='grey', border=F)
            }
            text(xlim[1]:xlim[2], rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]), '')[[1]], adj=c(0.5,1.3), cex=cex)
        }
        rect(starts - 0.5, 0, starts+0.5, h.del\$counts, col=col.del, border=h.del\$counts>0)
    } else if (force.plots) {
        plot(NA,xlim=c(-1,1), ylim=c(-1,1), axes=F, xlab=NA, ylab=NA, main=main)
        text(0, 0, labels=c('No deletions'))
    }

    ### Plot histogram for insertions only (if any)
    main=paste0('Location of the insertions (', label, ')')
    if (num.ins) {
        plot(h.ins\$counts, xlim=xlim, type='n', xlab=xlab.ins, ylab=ylab, main=main)
        if (nchar(ref.seq) > 0) {
            if (guide.start > 0 && guide.length > 0) {
                rect(starts[1] + guide.start - xlim[1] - 0.5, par('usr')[3]*0.975, starts[1] + guide.start - xlim[1] + guide.length - 0.5, 0, col='grey', border=F)
            }
            cex = 40/max(40, xlim[2]-xlim[1]+1)
            text(xlim[1]:xlim[2], rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]), '')[[1]], adj=c(0.5,1.3), cex=cex)
        }
        rect(starts.ins - 0.5, 0, starts.ins+0.5, h.ins\$counts, col=col.ins, border=h.ins\$counts>0)
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
    usr = par('usr');
    if (nchar(ref.seq) > 0) {
        factor = 10/100*40/max(40, xlim[2]-xlim[1]+1)
        par('usr' = c(usr[1], usr[2], usr[3]*(1+factor), usr[4]))
        if (guide.start > 0 && guide.length > 0) {
            rect(starts[1] + guide.start - xlim[1] - 0.5, -par('usr')[4]*factor, starts[1] + guide.start - xlim[1] + guide.length - 0.5, 0, col='grey', border=F)
        }
        ax.ticks = axTicks(2);
        axis(2, col=dark.col.del, col.axis=dark.col.del, at=subset(ax.ticks, ax.ticks>=0))
        rect(starts - 0.5, 0, starts+0.5, h.del\$counts, col=col.del, border=h.del\$counts>0)
        cex = 40/max(40, xlim[2]-xlim[1]+1)
        text(xlim[1]:xlim[2], rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]), '')[[1]], adj=c(0.5,1.3), cex=cex)
        par('usr' = c(usr[1], usr[2], usr[4], usr[3] *(1+factor)))
        axis(2, col=dark.col.ins, col.axis=dark.col.ins, at=subset(ax.ticks, ax.ticks>=0))
        rect(starts.ins - 0.5, 0, starts.ins+0.5, h.ins\$counts, col=col.ins, border=h.ins\$counts>0)
    } else {
        ax.ticks = axTicks(2);
        axis(2, col=dark.col.del, col.axis=dark.col.del, at=subset(ax.ticks, ax.ticks>=0))
        rect(starts - 0.5, 0, starts+0.5, h.del\$counts, col=col.del, border=h.del\$counts>0)
        par('usr' = c(usr[1], usr[2], usr[4], usr[3]))
        axis(2, col=dark.col.ins, col.axis=dark.col.ins, at=subset(ax.ticks, ax.ticks>=0))
        rect(starts.ins - 0.5, 0, starts.ins+0.5, h.ins\$counts, col=col.ins, border=h.ins\$counts>0)
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
    # Add ref sequence if available
    if (nchar(ref.seq) > 0) {
        if (guide.start > 0 && guide.length > 0) {
            rect(guide.start, par('usr')[3]*0.975, guide.start + guide.length, 0, col='grey', border=F)
        }
        cex = 40/max(40, xlim[2]-xlim[1]+1)
        text(xlim[1]:(xlim[2]-1)+0.5, rep(0, xlim[2]-xlim[1]+1), strsplit(substr(ref.seq, xlim[1], xlim[2]-1), '')[[1]], adj=c(0.5,1.3), cex=cex)
    }
    # add light blue rectangles with white background
    rect(min(range\$from):max(range\$to), 0, min(range\$from):max(range\$to)+1,
        del[min(range\$from):max(range\$to)], col=col.del, border='white')
    # add a black line around the profile of deletion frequencies
    lines(c(min(range\$from), min(range\$from):(max(range\$to)+1),min(range\$from)),
        c(0, del[min(range\$from):max(range\$to)], 0, 0), type='s', col='black')
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
        mtext('Complex: Deletion and insertion in the same sequence (or other mutations)', side=1, line=2)
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
        mtext(paste0('Complex: n = ', num.com, '; ', perc.com), side=1, line=2, adj=0, at=0.1)
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
#  FUNCTION plot.fastqc
# =============================================================================
#  This method plots a preview of the FastQC analysis
# =============================================================================
plot.fastqc <- function(fastqc.dir, title='FastQC') {
    library(png);
    mar = par('mar')
    par('mar' = c(0.2, 0.2, mar[3], 0.2))

    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, axes = F, main = title);

    summary <- read.table(paste0(fastqc.dir, '/summary.txt'), sep = '\t')
    tick <- readPNG(paste0(fastqc.dir, '/Icons/tick.png'))
    warning <- readPNG(paste0(fastqc.dir, '/Icons/warning.png'))
    error <- readPNG(paste0(fastqc.dir, '/Icons/error.png'))

    plot.flag <- function(test, x, y) {
        if (summary[summary\$V2 == test, 1] == 'PASS') {
            rasterImage(tick, x, y + 0.01, x + 0.04, y + 0.05)
        } else if (summary[summary\$V2 == test, 1] == 'WARN') {
            rasterImage(warning, x, y + 0.01, x + 0.04, y + 0.05)
        } else if (summary[summary\$V2 == test, 1] == 'FAIL') {
            rasterImage(error, x, y + 0.01, x + 0.04, y + 0.05)
        }
        text(x + 0.03, y + 0.025, test, pos = 4, cex = 0.75)
    }

    test <- 'Basic Statistics';
    plot.flag(test, 0, 0.95)
    data <- readLines(paste0(fastqc.dir, '/fastqc_data.txt'))
    line = 0
    lines(c(0, 0.3), c(0.915, 0.915))
    for (i in (which(data == '#Measure\tValue')[1] + 1):(which(data == '>>END_MODULE')[1] - 1)) {
        labels <- strsplit(data[i], '\t')[[1]]
        text(0, 0.9 - line, labels[1], pos = 4, cex = 0.35)
        text(0.16, 0.9 - line, labels[2], pos = 4, cex = 0.35)
        line <- line + 0.03
        lines(c(0, 0.3), c(0.915 - line, 0.915 - line))
    }
    lines(c(0, 0), c(0.915, 0.915 - line))
    lines(c(0.16, 0.16), c(0.915, 0.915 - line))
    lines(c(0.3, 0.3), c(0.915, 0.915 - line))

    test <- 'Per base sequence quality';
    plot.flag(test, 0.35, 0.95)
    png <- readPNG(paste0(fastqc.dir, '/Images/per_base_quality.png'))
    rasterImage(png, 0.35, 0.6, 0.65, 0.95)

    test <- 'Per tile sequence quality';
    plot.flag(test, 0.7, 0.95)
    png <- readPNG(paste0(fastqc.dir, '/Images/per_tile_quality.png'))
    rasterImage(png, 0.7, 0.6, 1, 0.95)

    test <- 'Per sequence quality scores';
    plot.flag(test, 0, 0.5)
    png <- readPNG(paste0(fastqc.dir, '/Images/per_sequence_quality.png'))
    rasterImage(png, 0,    0.15,  0.3, 0.5)

    test <- 'Per base N content';
    plot.flag(test, 0.35, 0.5)
    png <- readPNG(paste0(fastqc.dir, '/Images/per_base_n_content.png'))
    rasterImage(png, 0.35, 0.15, 0.65, 0.5)

    test <- 'Adapter Content';
    plot.flag(test, 0.7, 0.5)
    png <- readPNG(paste0(fastqc.dir, '/Images/adapter_content.png'))
    rasterImage(png, 0.7,  0.15,    1, 0.5)

    text(0.5, 0.05, paste0('See ', fastqc.dir, '.html for more details'), cex=0.6)

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
";

    my ($name, $path, $suffix) = fileparse($data_file);
    my $fastqc_R1_file = $data_file;
    $fastqc_R1_file =~ s/bam.data.txt/R1_val_1_fastqc/;
    if (-e "$fastqc_R1_file.zip") {
        system("unzip -q -u $fastqc_R1_file.zip -d $path");
        print R "

    plot.fastqc('$fastqc_R1_file', 'FastQC - R1');
";
    }

    my $fastqc_R2_file = $data_file;
    $fastqc_R2_file =~ s/bam.data.txt/R2_val_2_fastqc/;
    if (-e "$fastqc_R2_file.zip") {
        system("unzip -q -u $fastqc_R2_file.zip -d $path");
        print R "

    plot.fastqc('$fastqc_R2_file', 'FastQC - R2');
";
    }

print R "
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


sub test_merge_reads {
    my ($read1, $read2, $merged_read);

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "80M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads all matches");

    $read1 = {
        'seq' =>  "AAATTACCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "50M",
        'md' => "MD:Z:5A44"};
    $read2 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTACCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "80M",
                'md' => "MD:Z:5A74",
            },
            '5prime' => {
                'seq' =>  "AAATTACCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:5A24",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads mismatch 5'");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTACCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "50M",
        'md' => "MD:Z:35A14"};
    $read2 = {
        'seq' =>  "AAATTACCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:5A44"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTACCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "80M",
                'md' => "MD:Z:35A44",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTACCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:5A14",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads mismatch overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTAACGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "50M",
        'md' => "MD:Z:35A0A13"};
    $read2 = {
        'seq' =>  "AAATTAACGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:5A0A43"};
    $merged_read = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTAACGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "80M",
        'md' => "MD:Z:35A0A43"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTAACGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "80M",
                'md' => "MD:Z:35A0A43",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTAACGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:5A0A13",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads 2bp mismatch in the overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "10M1I40M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $merged_read = {
        'seq' =>  "AAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "10M1I70M",
        'md' => "MD:Z:80"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "10M1I70M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGTAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "10M1I20M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads insertion 5'");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "40M1I10M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "AAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "10M1I40M",
        'md' => "MD:Z:50"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "40M1I40M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGTAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "10M1I10M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads insertion 5'");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M1I20M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "30M1I50M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGT",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M1I",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads insertion just before overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M1I20M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "TAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "1I50M",
        'md' => "MD:Z:50"};
    $merged_read = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "30M1I50M",
        'md' => "MD:Z:80"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "30M1I50M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "TAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "1I20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads insertion as 1st bp overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTTAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M3I20M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "TTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "2I50M",
        'md' => "MD:Z:50"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "30M3I50M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGT",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M1I",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "TTAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "2I20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads partial insertion at start of overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTT",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M3I",
        'md' => "MD:Z:30"};
    $read2 = {
        'seq' =>  "TTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "2I50M",
        'md' => "MD:Z:50"};
    $merged_read = {
        'error' =>  $STAT_NO_OVERLAP};
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads overlap on insertion only");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTT",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M2I",
        'md' => "MD:Z:30"};
    $read2 = {
        'seq' =>  "CGGTTTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 28,
        'cigar' => "3M3I50M",
        'md' => "MD:Z:53"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "30M3I50M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCC",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 27,
                'cigar' => "27M",
                'md' => "MD:Z:27",
            },
            'overlap' => {
                'seq' =>  "CGGTT",
                'qual' => "IIIII",
                'start' => 28,
                'end' => 30,
                'cigar' => "3M2I",
                'md' => "MD:Z:3",
            },
            '3prime' => {
                'seq' =>  "TAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 80,
                'cigar' => "1I50M",
                'md' => "MD:Z:50",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads partial insertion at end of overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTT",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "30M3I",
        'md' => "MD:Z:30"};
    $read2 = {
        'seq' =>  "CGGTTTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 28,
        'cigar' => "3M3I50M",
        'md' => "MD:Z:53"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGTTTAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "30M3I50M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCC",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 27,
                'cigar' => "27M",
                'md' => "MD:Z:27",
            },
            'overlap' => {
                'seq' =>  "CGGTTT",
                'qual' => "IIIIII",
                'start' => 28,
                'end' => 30,
                'cigar' => "3M3I",
                'md' => "MD:Z:3",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 80,
                'cigar' => "50M",
                'md' => "MD:Z:50",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads partial insertion at end of overlap");


    $read1 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGATAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "31M1I19M",
        'md' => "MD:Z:50"};
    $read2 = {
        'seq' =>  "ATAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "1M1I49M",
        'md' => "MD:Z:50"};
    $merged_read = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGATAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "31M1I49M",
        'md' => "MD:Z:80"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGATAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "31M1I49M",
                'md' => "MD:Z:80",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
            'overlap' => {
                'seq' =>  "ATAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "1M1I19M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads insertion in the overlap");

    $read1 = {
        'seq' =>  "AAATTCCCGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "9M1D40M",
        'md' => "MD:Z:9^G40"};
    $read2 = {
        'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 31,
        'cigar' => "50M",
        'md' => "MD:Z:50"};
    $merged_read = {
        'seq' =>  "AAATTCCCGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "9M1D70M",
        'md' => "MD:Z:9^G70"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "9M1D70M",
                'md' => "MD:Z:9^G70",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 30,
                'cigar' => "9M1D20M",
                'md' => "MD:Z:9^G20",
            },
            'overlap' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIII",
                'start' => 31,
                'end' => 50,
                'cigar' => "20M",
                'md' => "MD:Z:20",
            },
            '3prime' => {
                'seq' =>  "AAATTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads deletion 5'");

    $read1 = {
        'seq' =>  "AAATTCCCGGBBBTTCCCGGDDDTTCCCTTCCCGGEEETTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'cigar' => "28M5D17M",
        'md' => "MD:Z:28^GGAAA17"};
    $read2 = {
        'seq' =>  "TTCCCGGEEETTCCCGGFFFTTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 34,
        'cigar' => "47M",
        'md' => "MD:Z:47"};
    $merged_read = {
        'seq' =>  "AAATTCCCGGBBBTTCCCGGDDDTTCCCTTCCCGGEEETTCCCGGFFFTTCCCGGAAATTCCCGGAAATTCCCGG",
        'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        'start' => 1,
        'end' => 80,
        'cigar' => "28M5D47M",
        'md' => "MD:Z:28^GGAAA47"};
    $merged_read = {
            'merged' => {
                'seq' =>  "AAATTCCCGGBBBTTCCCGGDDDTTCCCTTCCCGGEEETTCCCGGFFFTTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 80,
                'cigar' => "28M5D47M",
                'md' => "MD:Z:28^GGAAA47",
            },
            '5prime' => {
                'seq' =>  "AAATTCCCGGBBBTTCCCGGDDDTTCCC",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 1,
                'end' => 33,
                'cigar' => "28M5D",
                'md' => "MD:Z:28^GGAAA",
            },
            'overlap' => {
                'seq' =>  "TTCCCGGEEETTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIII",
                'start' => 34,
                'end' => 50,
                'cigar' => "17M",
                'md' => "MD:Z:17",
            },
            '3prime' => {
                'seq' =>  "FFFTTCCCGGAAATTCCCGGAAATTCCCGG",
                'qual' => "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
                'start' => 51,
                'end' => 80,
                'cigar' => "30M",
                'md' => "MD:Z:30",
            },
        };
    is_deeply(merge_reads($read1, $read2), $merged_read, "merge_reads deletion before overlap");

}

