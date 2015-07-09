#! /usr/bin/env perl
use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use File::Basename;
use IPC::Cmd qw[run];


=pod

=head1 NAME

crispr-parsr - Software to parse and analyse deletions (and insertions) from a CRISPR resequencing experiment

=head1 SYNOPSIS

crispr-parsr.pl [--wt ref_seq.fa] [--samples merge.txt] --input INPUT_DIR --output OUTPUT_DIR

=head1 DESCRIPTION

Briefly, the pipeline does the following:

=over 8

=item * Index the wild-type sequence for Bowtie2

=item * Merge FASTQ files according to the information in the samples file

=item * Trim and QC files with Trim Galore! (which uses cutadapt and FastQC internally)

=item * Align the reads with Bowtie2

=item * Parses the alignment to extract the deletions (and insertions) found in each sample

=back

=head1 OPTIONS

Please note that options names may be abbreviated to uniqueness, case does not matter, and a
single dash is sufficient

=over 8

=item B<-h>|B<--help>

Prints the full help.

=item B<--input INPUT_DIR>

The directory where all the FASTQ files are. By default, crispr-parsr will look for a FASTA file
(for the wild-type sequence) and a TXT file (for the samples definition) in this directory as well.

=item B<--output OUTPUT_DIR>

All files with intermediary data and final plots will be created in this directory. There will be a
README.txt file explaining the content of each file.

=item B<--ref-seq>|B<--wt-seq> ref_seq.fa>

This is a simple FASTA file with one sequence only (typically a short one) corresponding to the expected
wild-type sequence, before editing has ocurred.

If the filename is not provided, crispr-parsr will look for a FASTA file within the INPUT_DIR. If either
none or more than one file with the .fa extension exists in the INPUT_DIR, this will fail.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

    Example:
    ---------------------------------------
    >sample_amplicon1_ref_123456
    AACAGTGGGTCTTCGGGAACCAAGCAGCAGCCATGGGAGGTCTGGCTGTGTTCAGGCT
    CTGCTCGTGTAGATTCACAGCGCGCTCTGAACCCCCGCTGAGCTACCGATGGAAGAGG
    AGGAGGTCCTACAGTCGGAGATTCACAGCGCGCTCTGAACCACTTTCAGGAGACTCGA
    CTATTATGACTTATACGCGATA
    ---------------------------------------

=item B<--guide-seq> guide_seq.fa>

This is a simple FASTA file (header is optional) with one sequence only (typically a very short one)
corresponding to the guide sequence. It must match perfectly the reference sequence.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

    Example:
    ---------------------------------------
    AGGAGGTCCTA
    ---------------------------------------

=item B<--samples merge.txt>

This file contains the relation of FASTQ files for each sample. As the pipeline has been designed for
paired end reads, you need to specify the list of R1 and R2 fastq files in consecutive lines. You can
have more than one set of PE reads for each sample. For this, you need to separate your filenames with
a [tab] (you can use Excel and save the file in tab-separated format). Please use the same order for
R1 and R2 files within a sample.

You can also specify a label for each sample (recommended).

If the filename is not provided, crispr-parsr will look for a TXT file within the INPUT_DIR. If either
none or more than one file with the .txt extension exists in the INPUT_DIR, this will fail.

You can either provide the full path to the file or simply the name of the file if it is located
in the INPUT_DIR.

    Example (with labels):
    ---------------------------------------
    Mock:
    12345A01_S12_L001_R1_001.fastq
    12345A01_S12_L001_R2_001.fastq
    Run1:
    12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
    12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
    Run2:
    12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
    12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
    ---------------------------------------

    Example (without labels):
    ---------------------------------------
    12345A01_S12_L001_R1_001.fastq
    12345A01_S12_L001_R2_001.fastq
    12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
    12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
    12345D01_S22_L001_R1_001.fastq  12345F02_S41_L001_R1_001.fastq
    12345D01_S22_L001_R2_001.fastq  12345F02_S41_L001_R1_001.fastq
    ---------------------------------------

=back

=head1 Requirements

=over 8

=item B<TrimGalore!>: http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

Successfully tested with v0.3.3 and v0.4.0

=item B<FastQC>: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Successfully tested with v0.11.2 and v.0.11.3

=item B<cutadapt>: https://code.google.com/p/cutadapt/

Successfully tested with v1.8.1

=item B<samtools>: http://www.htslib.org

Successfully tested with v1.2

=item B<Bowtie2>: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Successfully tested with v2.2.3

=item B<R>: http://www.r-project.org

Successfully tested with v3.1.1

=back

Make sure that these softwares are available on the current path when running this tool.
Alternatively, you can hard-code the path at the top of the crispr-parsr.pl script. However
this does not work for fastqc and cutadapt which are called from within TrimGalore!.

=head1 INTERNAL METHODS

The rest of the documentation refers to the internal methods within this software and is
intended for developers only.

=cut


my $VERSION = "0.1a";

my $help;
my $verbose;
my $input_dir;
my $wt_seq_file;
my $guide_seq_file;
my $merge_file;
my $output_dir;
my $trim_galore = "trim_galore";
my $samtools = "samtools";
my $bowtie2_build_exe = "bowtie2-build";
my $bowtie2_exe = "bowtie2";
my $plotr = dirname(__FILE__)."/crispr-plotr.pl";

GetOptions(
    "help" => \$help,
    "verbose" => \$verbose,
    "input=s" => \$input_dir,
    "wt-seq|ref-seq|wt_file=s" => \$wt_seq_file,
    "guide_seq|guide-seq=s" => \$guide_seq_file,
    "samples|merge_file=s" => \$merge_file,
    "output=s" => \$output_dir,    
);

if ($help) {
    pod2usage(-verbose=>2);
}

if (!$input_dir or !$output_dir) {
    pod2usage(-verbose=>1);
}


if (!-d $output_dir) {
    mkdir($output_dir) or die "Cannot create output directory ($output_dir): $!\n";
}

create_readme_file($output_dir);

my $this_wt_seq_file = get_file_from_input_dir($input_dir, ".fa", $wt_seq_file);

my $this_guide_seq_file = get_file_from_input_dir($input_dir, ".fa", $guide_seq_file) if ($guide_seq_file);

my $bowtie2_index_file = index_genome($bowtie2_build_exe, $this_wt_seq_file, $output_dir);

my $this_merge_file = get_file_from_input_dir($input_dir, ".txt", $merge_file);

my $merged_files = merge_fastq_files($input_dir, $this_merge_file, $output_dir);

my $validated_files = run_trim_galore($trim_galore, $merged_files, $output_dir);

my $bam_files = run_bowtie2($bowtie2_exe, $bowtie2_index_file, $validated_files, $output_dir);

run_plotr($plotr, $bam_files, $this_wt_seq_file, $this_guide_seq_file, $output_dir);


=head2 create_readme_file

  Arg[1]        : string $output_dirname
  Example       : create_readme_file($output_dirname);
  Description   : Writes a README.txt file on the output directory. The content is mostly hard-coded
  Returns       : N/A
  Exceptions    : Dies if fails to create file

=cut

sub create_readme_file {
    my ($output_dir) = @_;

    open(README, ">$output_dir/README.txt") or die "Cannot open README.txt file in output directory ($output_dir)\n";

    print README qq{This directory contains the result of CRISPR-PARSR v$VERSION from the input in <$input_dir>


Briefly, the pipeline does the following:
 1. Index the wild-type sequence for using Bowtie2
 2. Merge FASTQ files according to the information in the merging file
 3. Trim and QC files with Trim Galore! (which uses cutadapt and FastQC internally)
 4. Align the reads with Bowtie2
 5. Parses the alignment to extract the deletions (and insertions) found in each sample

The information files (see below) contain the actual command line executed in each case.


Here is a brief description of the files you will find in here:

* Reference files

 - bt2/                             Directory containing the indexed wild-type sequence for Bowtie2


* Intermediary files

 - "sample".R1.fastq                Reads for this sample (first pair, after merging)
 - "sample".R2.fastq                Reads for this sample (second pair, after merging)

 - "sample".R1.fastq_trimming_report.txt    Trimming report for reads in this sample (first pair)
 - "sample".R2.fastq_trimming_report.txt    Trimming report for reads in this sample (second pair)

 - "sample".R1_val_1.fq             Validated reads for this sample (first pair)
 - "sample".R2_val_2.fq             Validated reads for this sample (second pair)

 - "sample".R1_val_1_fastqc.html    Web page with FastQC output for reads in this sample (first pair)
 - "sample".R1_val_1_fastqc.zip     ZIP file with FastQC output for reads in this sample (first pair)
 - "sample".R2_val_2_fastqc.html    Web page with FastQC output for reads in this sample (second pair)
 - "sample".R2_val_2_fastqc.zip     ZIP file with FastQC output for reads in this sample (second pair)

 - "sample".bam                     Alignments for this sample (both pairs)

 - "sample".bam.data.txt            Table with deletions in this sample


* Information files (actual command line + verbose output of these processes)

 - bowtie2-build.out.txt            Info from the bowtie2-build process (indexing the wild-type seq)
 - trim_galore."sample".out.txt     Info from Trim Galore! for this sampe (trimming and QC)
 - bowtie2."sample".out.txt         Info from the bowtie2 process for this sample (mapping the reads)
 - plotr."sample".out.txt           Info from the parsing and plotting process for this sample (extracting insertion and deletions)


* Results

 - "sample".out.txt                 Stats on WT/insert/deletion/filtered reads + most common deletions
 - "sample".out.pdf                 Plots showing the location and size of the deletions (PDF format)
 - "sample".out.XX.png              Plots showing the location and size of the deletions (PNG format)
 - "sample".out.XX.svg              Plots showing the location and size of the deletions (SVG format)
};
    close(README);

}


=head2 get_file_from_input_dir

  Arg[1]        : string $input_dir
  Arg[2]        : string $extension
  Arg[3]        : string $filename (optional)
  Example       : my $file = get_file_from_input_dir(".", ".txt");
  Example       : my $file = get_file_from_input_dir(".", ".txt", "file.txt");
  Description   : Search for a file with the provided extension in the input directory. This is
                  typically meant for cases where just one file of that type (extension) is found
                  in the input directory. The 3rd argument can be used in other cases to specify
                  which file to look for. In that case, the method confirms that the file exists.
  Returns       : string filename with full/relative path
  Exceptions    : Dies if no particular filename has been specified and more than 1 file found in
                  the input directory. Dies if no file is found.

=cut

sub get_file_from_input_dir {
    my ($input_dir, $extension, $file) = @_;
    my $this_file;

    if (!$input_dir) {
        die "Cannot read \*$extension file from unknown directory\n";
    }

    if ($file) {

        # In case the file has been provided
        if (-e $file) {
            $this_file = $file;
        } elsif (!-e "$input_dir/$file") {
            die "Cannot find \*$extension file $file in the input directory ($input_dir)\n";
        } else {
            $this_file = "$input_dir/$file";
        }

    } else {

        # Otherwise, try to find the (hopefully) unique file with that extension in the input directory
        opendir(INPUT, $input_dir) or die;
        my @files = grep {/$extension$/} readdir(INPUT);
        closedir(INPUT);

        if (!@files) {
            die "Cannot find any \*$extension file in the input directory ($input_dir)\n";
        }
        if (@files > 1) {
            die "More than one \*$extension file in the input directory ($input_dir)\n";
        }
        $this_file = "$input_dir/".$files[0];

    }
    
    return $this_file;
}


=head2 index_genome

  Arg[1]        : string $bowtie_build (the name of the executable, incl. the full path if needed)
  Arg[2]        : string $fasta_filename
  Arg[3]        : string $output_dirname
  Example       : my $bt2_index = index_genome($bowtie_build, $fasta_filename, $output_dirname);
  Description   : Index a fasta file for bowtie2
  Returns       : (string) Index filename prefix for its use with bowtie2
  Exceptions    : 

=cut

sub index_genome {
    my ($bowtie2_build_exe, $this_wt_seq_file, $output_dir) = @_;

    mkdir("$output_dir/bt2") unless (-d "$output_dir/bt2");

    my $command = [$bowtie2_build_exe, $this_wt_seq_file, "$output_dir/bt2/wt"];
    my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
    if (!$ok) {
        die "ERROR while running bowtie2-build: $err\n";
    }
    # Saves the command line and both stdout and stderr buffers (i.e. full buffer).
    open(OUT, ">$output_dir/bowtie2-build.out.txt") or die;
    print OUT "CMD: ", join(" ", @$command), "\n\n";
    print OUT "OUTPUT:\n";
    print OUT join("\n", @$full_buff);
    close(OUT);

    return "$output_dir/bt2/wt";
}


=head2 merge_fastq_files

  Arg[1]        : string $input_dir
  Arg[2]        : string $merge_filename
  Arg[3]        : string $output_dirname
  Example       : my $merged_files = merge_files($input_dir, $merge_filename, $output_dirname);
  Description   : Takes the list of files to be merged from the $merge_filename and merges them. The
                  name of the new files are returned. Note that this is intended for PE reads only.
  Returns       : hashref or arrays. ($label => [$merged_R1, $merged_R2])
  Exceptions    : Dies if problems with the input or the output.

=cut

sub merge_fastq_files {
    my ($input_dir, $merge_file, $output_dir) = @_;
    my $merged_files;
    my $merge_counter = 1;

    # Build a hash ($all_files) to quickly check if a file exists in the input dir. As a bonus, also
    # store the path to that file in the values of the hash.    
    my $all_files;
    opendir(INPUT, $input_dir) or die;
    map {$all_files->{$_} = "$input_dir/$_"} readdir(INPUT);
    closedir(INPUT);

    # Read UNIX, Windows or Mac file by slurping everything into a single string, substituting and splitting
    open(MERGE, $merge_file) or die;
    my $whole_text = join("\n", <MERGE>);
    $whole_text =~ s/[\r\n]+/\n/g;
    my @lines = split("\n", $whole_text);

    for (my $i = 0; $i < @lines - 1; $i += 2) {
        my $label;
        if ($lines[$i] =~ /^(\w+)\:$/ and $i < @lines - 2) {
            $label = $1;
            $i++;
        } else {
            $label = "sample.$merge_counter";
        }
        my $files1 = get_filenames_in_line($all_files, $lines[$i]);
        my $files2 = get_filenames_in_line($all_files, $lines[$i+1]);
        if (@$files1 != @$files2) {
            die "The number of FASTQ files on lines $i and ", ($i+1), " of the merge file ($merge_file) do not match\n";
        }
        $merge_counter++;
        my $merged_filename1 = "$output_dir/$label.R1.fastq";
        my $merged_filename2 = "$output_dir/$label.R2.fastq";
        my $command = ["cat", @$files1, ">", $merged_filename1];
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running cat: $err\n";
        }
        if (!-e $merged_filename1) {
            die "Cannot create merged file ($merged_filename1) in the output directory ($output_dir)\n";
        }
        $command = ["cat", @$files2, ">", $merged_filename2];
        ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running cat: $err\n";
        }
        if (!-e $merged_filename2) {
            die "Cannot create merged file ($merged_filename2) in the output directory ($output_dir)\n";
        }
        $merged_files->{$label} = [$merged_filename1, $merged_filename2];
    }

    return $merged_files;
}


=head2 run_trim_galore

  Arg[1]        : string $trim_galore (the name of the executable, incl. the full path if needed)
  Arg[2]        : hashref $merged_files ($label => [$fastq_R1, $fastq_R2])
  Arg[3]        : string $output_dirname
  Example       : my $validated_files = run_trim_galore($trim_galore, $merged_files, $output_dirname);
  Description   : Runs Trim Galore! on paired files and returns the name of the validated files
  Returns       : hashref or arrays. ($label => [$val_fastq_R1, $val_fastq_R2])
  Exceptions    : Dies if problems with the input or the output.
  
  TODO          : Check the output of FASTQC

=cut

sub run_trim_galore {
    my ($trim_galore, $merged_files, $output_dir) = @_;
    my $validated_files;

    while (my ($label, $this_pair_of_merged_files) = each %$merged_files) {
        die "The pair of merged files is not a pair\n" if (@$this_pair_of_merged_files != 2);
        my $command = [$trim_galore, "--fastqc", "-o", $output_dir, "--paired", @$this_pair_of_merged_files];
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running Trim Galore!: $err\n";
        }
        # Saves the command line and both stdout and stderr buffers (i.e. full buffer).
        # Note that it is not worth separating the STDOUT and STDERR buffers
        # in this case as Trim Galore! mixes them up.
        open(OUT, ">$output_dir/trim_galore.$label.out.txt") or die;
        print OUT "CMD: ", join(" ", @$command), "\n\n";
        print OUT "OUTPUT:\n";
        print OUT join("\n", @$full_buff);
        close(OUT);
        
        # Checks that the output files exist and store them in the array
        my ($validated_file1, $validated_file2) = @$this_pair_of_merged_files;
        $validated_file1 =~ s/.(fq|fastq|fa.gz|fastq.gz)$/_val_1.fq/;
        $validated_file2 =~ s/.(fq|fastq|fa.gz|fastq.gz)$/_val_2.fq/;
        if (!-e $validated_file1) {
            die "ERROR: Cannot find validated file $validated_file1\n";
        }
        if (!-e $validated_file2) {
            die "ERROR: Cannot find validated file $validated_file2\n";
        }

        my ($fastqc_zip_file1, $fastqc_zip_file2) = @$this_pair_of_merged_files;
        $fastqc_zip_file1 =~ s/.(fq|fastq|fa.gz|fastq.gz)$/_val_1_fastqc.zip/;
        $fastqc_zip_file2 =~ s/.(fq|fastq|fa.gz|fastq.gz)$/_val_2_fastqc.zip/;
        
        if (check_fastqc_zip_file($fastqc_zip_file1) and check_fastqc_zip_file($fastqc_zip_file2)) {
            $validated_files->{$label} = [$validated_file1, $validated_file2];
        }
    }
    
    return $validated_files;
}


=head2 check_fastqc_zip_file

  Arg[1]        : string $fastqc_zip_file (the name of the file, incl. the full path if needed)
  Example       : my $ok = check_fastqc_zip_file($fastqc_zip_file);
  Description   : Checks the results of FastQC stored in the summary.txt file. Some warnings and
                  failures are expected for this analysis, but other should be considered as errors.
  Returns       : boolean
  Exceptions    : Dies if cannot extract summary.txt from the file

=cut

sub check_fastqc_zip_file {
    my ($fastqc_zip_file) = @_;

    my $command = ["unzip", "-c", $fastqc_zip_file, "\*/summary.txt"];
    my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
    if (!$ok) {
        die "ERROR while unzipping $fastqc_zip_file!: $err\n";
    }
    
    my @summary_lines = split(/[\r\n]/, join("\n", @$stdout_buff));

    foreach my $this_line (@summary_lines) {
        if ($this_line =~ /^(PASS|WARN|FAIL)/) {
            my ($flag, $test, $file) = split("\t", $this_line);
            if ((($test eq "Basic Statistics") or
                ($test eq "Per base sequence quality") or
                ($test eq "Per tile sequence quality") or
                ($test eq "Per sequence quality scores") or
                ($test eq "Per base N content") or
                ($test eq "Adapter Content")) and $flag ne "PASS") {
                    print STDERR "FastQC failure: $file has a <$flag> for test <$test>\n";
                    print STDERR " -- The processing for this sample will stop here.\n";
                    return 0;
            }
        }
    }

    return "true";
}    

=head2 run_bowtie2

  Arg[1]        : string $bowtie2 (the name of the executable, incl. the full path if needed)
  Arg[2]        : string $bowtie2_index_filename_prefix (to be used in bowtie2 command line)
  Arg[3]        : hashref $fastq_files ($label => [$fastq_R1, $fastq_R2])
  Arg[4]        : string $output_dirname
  Example       : my $bam_files = run_bowtie2($bowtie2_exe, $bowtie2_index_file, $validated_files, $output_dir);
  Description   : Runs Bowtie2 on the input FASTQ files using the given index. Output is BAM-transformed
                  and stored in the output directory. The method returns the list of BAM files.
  Returns       : hashref of BAM filenames ($label => $bam_file)
  Exceptions    : Wrong format for the matrix or errors running Bowtie2

=cut

sub run_bowtie2 {
    my ($bowtie2_exe, $bowtie2_index_file, $validated_files, $output_dir) = @_;
    my $bam_files;

    while (my ($label, $this_pair_of_validated_files) = each %$validated_files) {
        my $this_bam_file = "$output_dir/$label.bam";
        if (@$this_pair_of_validated_files != 2) {
            die "Unexpected set of validated files: ".join(", ", @$this_pair_of_validated_files);
        }
        my ($file1, $file2) = @$this_pair_of_validated_files;
        my $command = [$bowtie2_exe, "-x", $bowtie2_index_file, "-1", $file1, "-2", $file2,
                        "|", $samtools, "view", "-Sb", "-",
                        ">", $this_bam_file];
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running bowtie2: $err\n";
        }
        # Saves the command line and the stderr buffer.
        # The stdout buffer is the SAM output captured in the command and stored in the bam file
        # using samtools.
        open(OUT, ">$output_dir/bowtie2.$label.out.txt") or die;
        print OUT "CMD: ", join(" ", @$command), "\n\n";
        print OUT "OUTPUT:\n";
        print OUT join("\n", @$stderr_buff);
        close(OUT);

        # Checks that the output file exists
        if (!-e $this_bam_file) {
            die "ERROR: Cannot find expected bam file $this_bam_file\n";
        }
        $bam_files->{$label} = $this_bam_file;
    }

    return $bam_files;
}


=head2 run_plotr

  Arg[1]        : string $plotr (the name of the Perl script with full path)
  Arg[2]        : hashref $bam_files ($label => $bam_file)
  Arg[3]        : string $wt_fasta_file
  Arg[4]        : string $output_dirname
  Example       : run_plotr($plotr, $bam_files, $wt_fasta_file, $output_dir);
  Description   : 
  Returns       : 
  Exceptions    : 

=cut

sub run_plotr {
    my ($plotr, $bam_files, $this_wt_seq_file, $this_guide_seq_file, $output_dir) = @_;

    while (my ($label, $this_bam_file) = each %$bam_files) {
        my $this_pdf_file = "$output_dir/$label.out.pdf";
        my $this_txt_file = "$output_dir/$label.out.txt";
        my $command = [$plotr, "--input", $this_bam_file, "--out", $this_pdf_file, "--ref_seq", $this_wt_seq_file,
                        "--label", $label];
        push(@$command, "--guide-seq", $this_guide_seq_file) if ($this_guide_seq_file);
        push(@$command, ">", $this_txt_file);
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running plotr: $err\n";
        }
        # Saves the command line and the stderr buffer.
        # The stdout buffer is the TXT output captured in the command.
        open(OUT, ">$output_dir/plotr.$label.out.txt") or die;
        print OUT "CMD: ", join(" ", @$command), "\n\n";
        print OUT "OUTPUT:\n";
        print OUT join("\n", @$stderr_buff);
        close(OUT);

        if (!-e $this_pdf_file) {
            die "ERROR: Cannot find expected PDF file $this_pdf_file\n";
        }
        if (!-e $this_txt_file) {
            die "ERROR: Cannot find expected TXT file $this_txt_file\n";
        }
    }
}


=head2 get_filenames_in_line

  Arg[1]        : hashref $all_files (key: filename; value: filename with full/relative path)
  Arg[2]        : string $line
  Example       : my $files = get_filenames_in_line($all_files, $line);
  Description   : Splits the line using any space character. Each string is a filename. Checks that
                  the filename exists in the $all_files hash and return the list of files with the
                  full/relative paths.
  Returns       : arrayref of filenames with full/relative paths
  Exceptions    : Dies if any filename is not found.
  Caller        : merge_files

=cut

sub get_filenames_in_line {
    my ($all_files, $line) = @_;
    my $files = [];

    $line =~ s/\s+/\t/;
    @$files = split("\t", $line);
    foreach my $this_file (@$files) {
        if (!$all_files->{$this_file}) {
            die "FASTQ file $this_file not found in input directory\n";
        }
        $this_file = $all_files->{$this_file};
    }

    return $files;
}


