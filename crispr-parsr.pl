#! /usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;
use IPC::Cmd qw[run];

my $VERSION = "0.1a";

my $help;
my $verbose;
my $input_dir;
my $wt_seq_file;
my $merge_file;
my $output_dir;
my $trim_galore = "/Users/javier/Downloads/trim_galore_v0-2/trim_galore";
my $bowtie2_build_exe = "bowtie2-build";
my $bowtie2_exe = "bowtie2";

my $desc = qq{crispr-parsr.pl --input INPUT_DIR --output OUTPUT_DIR
};

GetOptions(
    "help" => \$help,
    "verbose" => \$verbose,
    "input=s" => \$input_dir,
    "wt_file=s" => \$wt_seq_file,
    "merge_file=s" => \$merge_file,
    "output=s" => \$output_dir,    
);

if ($help) {
    print $desc;
    exit(0);
}

mkdir($output_dir) or die "Cannot create output directory ($output_dir): $!\n";

create_readme_file($output_dir);

my $this_wt_seq_file = get_file_from_input_dir($input_dir, ".fa", $wt_seq_file);

my $bowtie2_index_file = index_genome($bowtie2_build_exe, $this_wt_seq_file, $output_dir);

my $this_merge_file = get_file_from_input_dir($input_dir, ".txt", $merge_file);

my $merged_files = merge_fastq_files($input_dir, $this_merge_file, $output_dir);

my $validated_files = run_trim_galore($trim_galore, $merged_files, $output_dir);

my $bam_files = run_bowtie2($bowtie2_exe, $bowtie2_index_file, $validated_files, $output_dir);

run_parser("parse_crispr_pe_bam_file.pl", $bam_files, $this_wt_seq_file, $output_dir);


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
 - parser."sample".out.txt          Info from the parser process for this sample (extracting insertion and deletions)


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
  Returns       : arrayref or arrays, i.e. ref to a matrix, with the name of the merged files. Each
                  row is a pair of files (for PE reads) to treat together.
  Exceptions    : Dies if problems with the input or the output.

=cut

sub merge_fastq_files {
    my ($input_dir, $merge_file, $output_dir) = @_;
    my $merged_files;
    my $merge_counter = 0;

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
            $label = "merged.$merge_counter";
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
  Arg[2]        : arrayref $merged_files (see merge_files)
  Arg[3]        : string $output_dirname
  Example       : my $validated_files = run_trim_galore($trim_galore, $merged_files, $output_dirname);
  Description   : Runs Trim Galore! on paired files and returns the name of the validated files
  Returns       : arrayref or arrays, i.e. ref to a matrix, with the name of the validated files. Each
                  row is a pair of files (for PE reads) to treat together.
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
        $validated_files->{$label} = [$validated_file1, $validated_file2];
    }
    
    #TODO: Check the FastQC output
    return $validated_files;
}


=head2 run_bowtie2

  Arg[1]        : string $bowtie2 (the name of the executable, incl. the full path if needed)
  Arg[2]        : string $bowtie2_index_filename_prefix (to be used in bowtie2 command line)
  Arg[3]        : arrayref $fastq_files (matrix with pairs of fastq files to align with Bowtie2
  Arg[4]        : string $output_dirname
  Example       : my $bam_files = run_bowtie2($bowtie2_exe, $bowtie2_index_file, $validated_files, $output_dir);
  Description   : Runs Bowtie2 on the input FASTQ files using the given index. Output is BAM-transformed
                  and stored in the output directory. The method returns the list of BAM files.
  Returns       : listref of BAM filenames
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
                        "|", "samtools", "view", "-Sb", "-",
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


=head2 run_parser

  Arg[1]        : string $parser (the name of the Perl script with full path)
  Arg[2]        : arrayref $bam_files (arrayref of bam files from Bowtie2)
  Arg[3]        : string $wt_fasta_file
  Arg[4]        : string $output_dirname
  Example       : run_parser($parser, $bam_files, $wt_fasta_file, $output_dir);
  Description   : 
  Returns       : 
  Exceptions    : 

=cut

sub run_parser {
    my ($parser, $bam_files, $this_wt_seq_file, $output_dir) = @_;

    while (my ($label, $this_bam_file) = each %$bam_files) {
        my $this_pdf_file = "$output_dir/$label.out.pdf";
        my $this_txt_file = "$output_dir/$label.out.txt";
        my $command = ["perl", $parser, "--input", $this_bam_file, "--out", $this_pdf_file, "--ref_seq", $this_wt_seq_file,
                        "--label", $label, ">", $this_txt_file];
        my ($ok, $err, $full_buff, $stdout_buff, $stderr_buff) = run(command => $command);
        if (!$ok) {
            die "ERROR while running parser: $err\n";
        }
        # Saves the command line and the stderr buffer.
        # The stdout buffer is the TXT output captured in the command.
        open(OUT, ">$output_dir/parser.$label.out.txt") or die;
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


