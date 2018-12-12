#!/usr/bin/perl
use List::Util qw[min max];
use Getopt::Long;
use Pod::Usage;
use strict;
use Cwd;
use IPC::System::Simple qw(system systemx capture capturex $EXITVAL run);

=pod 

=head1 NAME
    
    
    Program: ChIP-Seq_analysis_pipeline.pl
    Authors: Alina Saiakhova (ars51@case.edu)
    Last Update: 07/30/2013

=head1 SYNOPSIS

     
 perl ChIP-Seq_analysis_pipeline.pl -f chip.fastq -d /full/path/to/genome/prefix [options]
  
    HELP
  					-h,--help		prints this message

  
    REQUIRED PARAMETERS
  					-f,--chip		ChIP-Seq treatment file (acceptable formats are: SRA, gzipped SRA, FASTQ, gzipped FASTQ, BAM, and SAM). 
  					-d ,--refdir		path to reference files  (i.e. /path/to/genomes/genome_file_prefix where genome_file_prefix.fai is a file; used to locate files by BWA, Bowtie, and SAMtools)
 
 
    OPTIONAL PARAMETERS 
	    General options:
					-o,--chipname 		output prefix for ChIP-Seq sample (default is ChIP-Seq_sample). Example: V410_H3K4me1_hg19_05232013
					-i,--input    		non-ChIP/input DNA file (acceptable formats are: SRA, gzipped SRA, FASTQ, gzipped FASTQ, BAM, and SAM)
  					-c,--controlname	output prefix for INPUT sample (default is control_sample). Example: V410_input_hg19_05232013

	    Adapter/Quality Trimming options:
					-t ,--type		sequencing technology: illumina or sanger (default is illumina). ASCII offsets of 33 and 64 will be used by FastX tools on Illumina and Sanger data, respectively
					-atrim		 	adapter trim data
  					--chipadapter		ChIP adapter sequence to trim (default is standard Illumina adapter: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT)	
					--inputadapter		Input adapter sequence to trim (default is standard Illumina adapter: AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT)					
					-qtrim 			quality filter data
					-q ,--qual 		phred quality score for quality filtering (default is 20)
					-l, --len		min read length to keep during  adapter trimming or quality filtering (default is 25)

	    Alignment/SAMtools options:
					--aligner		alignment program to use for aligning FASTQ data, bwa,  bowtie or bowtie2 (default is bowtie2)
					-k			report a max of <k> alignments (works for bowtie and bowtie2 only)						
					-p, -thread		number of threads to use for BWA/Bowtie/Bowtie2 alignments (default is 4)
					--normdup		Do not remove potential PCR duplicates prior to peak calling. This will skip the SAMtools PCR duplicates removal step and run MACS with --keep-dup=TRUE. 

					

	    Peak Calling, WIG, SGR options:
					--nomodel		skip peak model building while running MACS
					--extsize		shift tags to the peak center by this distance in bp, if --nomodel is selected; will be ignored otherwise ( default is 200)
					--genomesize 		effective genome size, default is hs (human genome); use MACS genome shortcuts or specify actual genome size
					--broad			option must be set when called peaks are expected to be broad
					--nowig			do not create wig tracks
					--chromsizes		chromosome sizes file (for bigWIG generation)
					
	    Evaluating ChIP-seq data:
					--frip			calculate fraction of reads in peaks (FRiP)

    
  Dependencies: FastX Toolkit, BWA/Bowtie/Bowtie2, SAMtools, macs14, fastq-dump
 
  Example: ./ChIP-Seq_pipeline.pl -f chip.fastq -i input.fastq -o ChIP_name -c Input_name --refdir /home1/genomes/hg18/hg18 -s 10 -p 4

=cut

use Getopt::Long;
use Pod::Usage;


my %options 		= (
	thread				=> 4,
    type    		=> 'illumina',
    qual    		=> 20,
    len     		=> 25,
    chipname     => 'ChIP-Seq_sample',
    controlname => 'Control_sample',
    chipadapter		=> 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT',
    inputadapter		=> 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTT',
    aligner => 'bowtie2',
    genomesize => 'hs',
    
);

GetOptions(\%options, 
    'help|h', 
    'chip|f=s', 
    'input|i=s', 
    'chipname|o=s',
    'controlname|c=s',
    'refdir|d=s',
    'thread|p=i',
    'type|t=s', 
    'qual|q=i', 
    'len|m=i',
    'inputadapter=s',
    'chipadapter=s',
    'aligner=s',
    'k=i',
    'atrim',
    'qtrim',
    'genomesize=s',
    'nowig',
    'nomodel',
    'extsize',
    'broad',
    'frip',
    'wignorma',
    'wignormc',
    'wignormmean',
    'chromsizes=s',
  'normdup', 
);

####DEFINE GLOBAL VARIABLES####
my $ALIGNER=$options{aligner};
my $CHIP_PREFIX=$options{chipname};
my $MACS_OUT_PREFIX=$CHIP_PREFIX."_".$ALIGNER;
my $INPUT_PREFIX=$options{controlname};
my $NUM_ALIGNED_CHIP_READS=0;
my $NUM_ALIGNED_INPUT_READS=0;




###END GLOBAL VARIABLES#######


######USAGE/HELP INFO############
my $log_file="$CHIP_PREFIX"."_log.txt";
unlink $log_file if -f $log_file;
pod2usage(-verbose=>99) if  $options{'help'};
stderr( "Usage: perl ChIP-Seq_analysis_pipeline.pl -f chip.fastq -d /full/path/to/genome/prefix [options]\n")  if ($options{chip} eq "" or $options{refdir} eq "" );
stderr ("Please specify one of the following aligners: bwa, bowtie2, or bowtie2\n") if($ALIGNER!~/^(bowtie|bowtie2|bwa)$/);
######END USAGE/HELP INFO##########





############################START MAIN PIPELINE###########################################
stdout( "-------------------ChIP-Seq Pipeline started (".(localtime).")-------------------\n\n");


## Step 1A: run FASTX tools, aligner, and SAMtools on TREATMENT file
my $chip_file_for_macs=$options{chip};
my $input_file_for_macs=$options{input};
$chip_file_for_macs=run_pre_macs_steps($options{chip}, $CHIP_PREFIX, $options{chipadapter}) if $chip_file_for_macs!~/.bam/;

## Step 1B (optional): run FASTX tools, aligner and SAMtools on CONTROL file, if applicable
$input_file_for_macs=run_pre_macs_steps($options{input}, $INPUT_PREFIX, $options{inputadapter}) if ($options{input} and $input_file_for_macs!~/.bam/);

##Step 2: run MACS
run_macs_step($chip_file_for_macs,$input_file_for_macs);

##Step 3 (optional): wig tracks normalization (deprecated for now)

##wignormA($chip_file_for_macs) if $options{wignorma} and $options{chromsizes};
##wignormC() if $options{wignormc} and $options{input} and $options{chromsizes};
##wignormMean() if $options{wignormmean} and $options{chromsizes};



##Step 4 (optional): quality assessment of ChIP-Seq data
frip_calc($MACS_OUT_PREFIX, $chip_file_for_macs) if $options{frip};

stdout(  "\n\n-------------------ChIP-Seq Pipeline finished (".(localtime).")-------------------\n");
############################END MAIN PIPELINE###########################################
sub stderr{   
    my $message=$_[0];
    eval {
        run("echo -n \"$message\" >> $log_file");
    };

    if ($@) {
        print "Something went wrong with writing STDERR to file - $@\n";
    }

    die $message;
}
sub stdout{
    my $message=$_[0];
    print $message;
       
eval {
        run("echo -n \"$message\" >> $log_file");
    };

    if ($@) {
        print "Something went wrong with writing STDOUT to file - $@\n";
    }

    
    
}

sub run_command {
    my $command=shift;
    my $error_message_if_fail=shift;
    chomp $error_message_if_fail;
    
     eval {
       my $output  = capture($command);
       stdout($output."\n");
    };

    if ($@) {
        stderr "$error_message_if_fail - $@\n";
    }
    
    
    
    
}
sub run_pre_macs_steps{
    
    my $file=$_[0];
    my $file_prefix=$_[1];
   my $adapter=$_[2];
    
    ##figure out file format and reformat, if necessary
    if($file=~m/(.+)\.gz|\.gzip/i){
	run_command("gunzip -c $file >> $1 2>> $log_file","Decompressing $file failed" );
	$file=$1;
    }
    if($file=~/\.sra/i){
	
	run_command("fastq-dump -Z $file >> $file_prefix.fastq 2>> $log_file", "fastq-dump on $file failed");
	$file="$file_prefix.fastq";
    }
    elsif($file=~/\.sam/i){
	run_command("samtools view -b -S  $file > $file_prefix.bam 2>> $log_file ", "conversion of $file to BAM file failed");
	$file="$file_prefix.bam";
    }

    else {
	#assumes it's either a FASTQ file or a BAM file
    }
    
    if ($file!~/\.bam|\.sam/){ ##if raw data file, then proceed with adapter/quality filtering and alignment
	
	my $q_shift;
	if($options{type} eq "illumina"){$q_shift=33; }elsif($options{type} eq "sanger"){$q_shift=64; }else{ stderr ("Unrecognized sequencing method: $options{type}\n");}

	$file=adapter_trim_file($file, $file_prefix, $q_shift, $adapter) if $options{atrim} ; ##ADAPTER TRIMMING (optional)
	$file=quality_filter_file($file, $file_prefix, $q_shift) if $options{qtrim} ; ##QUALITY FILTERING (optional)
	$file=align_sample($file, $file_prefix, $q_shift); ##ALIGN DATA
    }
    
    $file=run_samtools_steps($file, $file_prefix); ##SORT, INDEX FILES; REMOVE PCR DUPLICATES (requires for all types of input files)

    return $file;
}

sub run_macs_step{
    my $chip_sample_bam=$_[0];
    my $input_sample_bam=$_[1];
    
    my $scaling=my $wiggen=my $modelgen=my $shiftsize=my $control_str=my $keep_dup=my $broad="";
    
    if ($input_sample_bam ne ""){
    $scaling=determine_scaling_param($chip_sample_bam, $input_sample_bam);
    $control_str="-c $input_sample_bam";
    }
   $wiggen="-B --trackline" if !$options{nowig};
   $modelgen="--nomodel" if $options{nomodel};
   $shiftsize="--extsize=$options{extsize}" if $options{nomodel} and $options{extsize};
   $keep_dup="--keep-dup=TRUE" if $options{normdup};
   $broad="--broad" if $options{broad};
   
    my $macs_call_peak_string="macs2 callpeak -t $chip_sample_bam --name=$MACS_OUT_PREFIX --format=BAM -g $options{genomesize}  $keep_dup $control_str $wiggen $broad $modelgen $shiftsize $scaling 2>&1";
        stdout (">>>>>>>>>Running MACS2 call peak step...\n");
    run_command($macs_call_peak_string, "MACS2 callpeak failed");
        stdout (">>>>>>>>>Creating lambda normalized bigWIGs...\n");
    my $macs_bdgcmp_string="macs2 bdgcmp -t ".$MACS_OUT_PREFIX."_treat_pileup.bdg -c ".$MACS_OUT_PREFIX."_control_lambda.bdg -m ppois --o-prefix ".$MACS_OUT_PREFIX."_control_lambda_normed";
    run_command($macs_bdgcmp_string, "MACS2 bdgcmp failed...\n") if !$options{nowig};
    run_command("sort -k1,1 -k2,2n ".$MACS_OUT_PREFIX."_control_lambda_normed_ppois.bdg > sorted_".$MACS_OUT_PREFIX."_control_lambda_normed_ppois.bdg","bedGraph sorting failed...\n");
    run_command("bedGraphToBigWig  sorted_".$MACS_OUT_PREFIX."_control_lambda_normed_ppois.bdg $options{chromsizes} $MACS_OUT_PREFIX.bw","bedGraphToBigWig step failed...\n");


    

}
sub determine_scaling_param{
    my $chip_sample_bam=$_[0];
    my $input_sample_bam=$_[1];
    my $scaling="";
    
    stdout (  ">>>>>>>>>Determining number of aligned reads in $chip_sample_bam...\n");
    my $num_aligned_chip_reads=`samtools view -c -F 4 $chip_sample_bam `;
    $NUM_ALIGNED_CHIP_READS=$num_aligned_chip_reads;
    chomp $num_aligned_chip_reads;
    stdout (  "Found $num_aligned_chip_reads reads\n");
    
   stdout (  ">>>>>>>>>Determining number of aligned reads in $input_sample_bam...\n");
    my $num_aligned_input_reads=`samtools view -c -F 4 $input_sample_bam `;
    $NUM_ALIGNED_INPUT_READS=$num_aligned_input_reads;
    chomp $num_aligned_input_reads;
    stdout (  "Found $num_aligned_input_reads reads\n");
	    
    if($num_aligned_chip_reads >= $num_aligned_input_reads){
	stdout (  "WARNING: $chip_sample_bam has ".sprintf("%.2f", $num_aligned_chip_reads/$num_aligned_input_reads, )."X more reads than $input_sample_bam. The input dataset will be upscaled.\n");
	$scaling="--to-large";
	 }else{
	
	stdout ( "WARNING: $input_sample_bam has ".sprintf("%.2f", $num_aligned_input_reads/$num_aligned_chip_reads, )."X more reads than $chip_sample_bam. The input dataset will be downscaled.\n");
    }
    
    return $scaling;
}
sub adapter_trim_file{
    my $file=$_[0];
    my $file_prefix=$_[1];
    my $q_shift=$_[2];
    my $adapter=$_[3];
    
 stdout(  ">>>>>>>>>Adapter trimming  $file ($options{type} format)...\n");
run_command ("fastx_clipper -a $adapter -l $options{len} -v -Q$q_shift -i $file -o clipped_$file_prefix.fastq 2>&1", "Adapter trimming on $file failed" );
return "clipped_$file_prefix.fastq";
}

sub quality_filter_file{
    my $file=$_[0];
    my $file_prefix=$_[1];
    my $q_shift=$_[2];
    
stdout ( ">>>>>>>>>Quality filtering  $file ($options{type} format)...\n");
run_command("fastq_quality_trimmer -t $options{qual} -l $options{len} -v -Q$q_shift -i $file -o trimmed_$file_prefix.fastq 2>&1","Quality trimming on $file failed");
return "trimmed_$file_prefix.fastq";


    
    
}


sub align_sample{
    my $file=$_[0];
    my $file_prefix=$_[1];
    my $q_shift=$_[2];

    my $aln_num_param="";
    $aln_num_param=" -k $options{k}" if $options{k};
    
    stdout (">>>>>>>>>Aligning $file_prefix ($file) with $ALIGNER...\n");
    if($ALIGNER eq "bwa"){
	
   run_command ("bwa aln -t $options{thread} $options{refdir} $file >> $file_prefix.$ALIGNER.sai 2>> $log_file","Error: sai file cannot be created for $file");
    
    stdout ( ">>>>>>>>>Generating SAM file...\n");
   run_command ("bwa samse $options{refdir} $file_prefix.$ALIGNER.sai $file >> $file_prefix.$ALIGNER.sam 2>> $log_file", "Error: sam file cannnot be created for $file_prefix.$ALIGNER.sai");
    }
    elsif($ALIGNER eq "bowtie"){
	
   run_command ("bowtie $options{refdir} -q -m 1 -p $options{thread} --phred$q_shift-quals -S $file >> $file_prefix.$ALIGNER.sam 2>> $log_file", "Error: Bowtie alignment failed on $file");

    }
    elsif($ALIGNER eq "bowtie2"){
	
    run_command ("bowtie2 $options{refdir} -p $options{thread}  $aln_num_param --phred$q_shift -U $file -S  $file_prefix.$ALIGNER.sam 2>&1", "Error: Bowtie2 alignment failed on $file");    

	
    }
    
    return "$file_prefix.$ALIGNER.sam";
}
sub run_samtools_steps{
    my $file=$_[0];
    my $file_prefix=$_[1];
    
    stdout ( ">>>>>>>>>Generating BAM file...\n");
    run_command ("samtools view -b -S  $file > $file_prefix.$ALIGNER.bam 2>> $log_file ", "BAM file cannot be created from $file");


    if (!$options{normdup}){
    stdout (  ">>>>>>>>>Sorting BAM file...\n");
    run_command ("samtools sort $file_prefix.$ALIGNER.bam $file_prefix.sorted.$ALIGNER 2>&1", "$file_prefix.$ALIGNER.bam cannot be sorted");

    stdout ( ">>>>>>>>>Indexing sorted BAM file for faster access...\n");
    run_command ("samtools index $file_prefix.sorted.$ALIGNER.bam 2>&1", "file_prefix.sorted.$ALIGNER.bam cannot be indexed");
    
    stdout (  ">>>>>>>>>Removing potential PCR duplicates..\n");
    run_command ("samtools rmdup -S $file_prefix.sorted.$ALIGNER.bam $file_prefix.rmdup.$ALIGNER.bam 2>&1", "Potential PCR duplicates cannot be removed from $file_prefix.sorted.$ALIGNER.bam");
    
         return "$file_prefix.rmdup.$ALIGNER.bam";

    } else {
	
		return "$file_prefix.$ALIGNER.bam";
	
    }
    
    
}





sub frip_calc{
    my $macs_prefix=$_[0];
    my $chip_bam_file=$_[1];
    my $peak_file=$macs_prefix."_peaks.bed";
    
    $NUM_ALIGNED_CHIP_READS=`samtools view -c -F 4 $chip_bam_file` if $NUM_ALIGNED_CHIP_READS==0;
    
    stdout (  ">>>>>>>>>Calculating FRiP based on $peak_file and $chip_bam_file...\n") ;
    my $frip_report=`cat $peak_file | awk '(\$2<0){\$2=0}{printf \"%s\\t%d\\t%d\\n\" , \$1, \$2, \$3}' | coverageBed -abam -counts -a $chip_bam_file -b - | cut -f 1-4 | awk '{total=total+\$4}END{print total, $NUM_ALIGNED_CHIP_READS,total/$NUM_ALIGNED_CHIP_READS}'`;
    chomp $frip_report;
    open FRIP_FILE, ">", "$macs_prefix"."_frip_report.txt";
    my @frip_report_arr=split(/\s+|\t/, $frip_report);
    printf FRIP_FILE ("There are %s reads in peaks out of a total of %s aligned reads. The FRIP score is %s\n", @frip_report_arr);
    
}

sub wignormMean{
stdout (  ">>>>>>>>>Normalizing treatment dataset $CHIP_PREFIX wiggle tracks by mean WIG signal...\n");
`gunzip $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*.gz `;
my $norm_by="mean";
my ($output_prefix, $chrom_sizes_file)=( $MACS_OUT_PREFIX."_".$norm_by."_normed", $options{chromsizes});
`cat $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*chr*.wig > pre_normed_$output_prefix.wig`;
my $file_to_norm="pre_normed_$output_prefix.wig";

open(IN, "<$file_to_norm") || die ("Can't open $file_to_norm\n");
my ($running_sum, $signal_line_count)=(0,0);
while(my $line=<IN>){
chomp $line;
if($line!~m/track|variable/){
my ($coord, $signal)=split(/\t|\s+/, $line);
$running_sum+=$signal;
$signal_line_count++;
}
}
my $norm_factor=$running_sum/$signal_line_count;
print "Normalizing $file_to_norm by $norm_by = $norm_factor...\n";
unlink($output_prefix.".wig");
system(" head -1 $file_to_norm> $output_prefix.wig");
run_command(" cat $file_to_norm | grep -v wiggle | awk ' {if (\$0 ~ /Step/){print \$0 }else {printf \"%d\t%0.3f\\n\", \$1,\$2/$norm_factor}}' >> $output_prefix.wig","norming failed");
run_command("wigToBigWig  $output_prefix.wig $chrom_sizes_file  $output_prefix.bw -clip 2>&1","wigToWig failed");
unlink("$output_prefix.wig");
unlink($file_to_norm);
`gzip $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*.wig `;

}
 
sub wignormC {
   stdout (  ">>>>>>>>>Normalizing treatment dataset $CHIP_PREFIX against $INPUT_PREFIX using wignorm...\n");
  `gunzip $MACS_OUT_PREFIX"_MACS_wiggle"/*/*.gz `;  
      my %genome_shortcuts_to_size=('hs'=>'2700000000','mm'=>'1870000000','ce'=>'90000000','dm'=>'120000000');
      my ($genome_for_wignorm, $output_prefix, $chrom_sizes_file)=($options{genomesize},$MACS_OUT_PREFIX."_wignorm", $options{chromsizes});
      $genome_for_wignorm=$genome_shortcuts_to_size{$options{genomesize}} if($genome_shortcuts_to_size{$options{genomesize}});
       my $whole_chip_wig="$CHIP_PREFIX.wig";
       my $whole_control_wig="$CHIP_PREFIX"."_control.wig";
     `cat $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*.wig > $whole_chip_wig`;
    `cat $MACS_OUT_PREFIX"_MACS_wiggle"/control/*.wig > $whole_control_wig`;
run_command("wignorm -t $whole_chip_wig -c $whole_control_wig --gsize=$genome_for_wignorm --name=$output_prefix 2>&1", "Wignorm failed");
run_command("wigToBigWig ".$output_prefix."_scores.wig  $chrom_sizes_file  $output_prefix.bw -clip");
`rm $whole_chip_wig $whole_control_wig $output_prefix"_scores.wig" `;
`gzip $MACS_OUT_PREFIX"_MACS_wiggle"/*/*.wig `;  

    
}

sub wignormA {
my $chip_bam_file=$_[0];
stdout (  ">>>>>>>>>Normalizing treatment dataset $CHIP_PREFIX wiggle tracks by the total number of aligned reads (in million)...\n");
`gunzip $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*.gz `;
my $num_aln_reads=`samtools view -c -F 4 $CHIP_PREFIX.rmdup.$ALIGNER.bam`;
my $num_aln_reads_in_mil=$num_aln_reads/1000000;
my $norm_by="aln";
my ($output_prefix,  $chrom_sizes_file)=( $MACS_OUT_PREFIX."_".$norm_by."_normed", $options{chromsizes});
`cat $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*chr*.wig > pre_normed_$output_prefix.wig`;
my $file_to_norm="pre_normed_$output_prefix.wig";

my $norm_factor=$num_aln_reads_in_mil;
stdout ("Normalizing $file_to_norm by $norm_by = $norm_factor...\n");
unlink($output_prefix.".wig");
system(" head -1 $file_to_norm> $output_prefix.wig");
run_command(" cat $file_to_norm | grep -v wiggle | awk ' {if (\$0 ~ /Step/){print \$0 }else {printf \"%d\t%0.3f\\n\", \$1,\$2/$norm_factor}}' >> $output_prefix.wig","norming failed");
run_command("wigToBigWig  $output_prefix.wig $chrom_sizes_file  $output_prefix.bw -clip 2>&1","wigToWig failed");
unlink("$output_prefix.wig");
unlink($file_to_norm);
`gzip $MACS_OUT_PREFIX"_MACS_wiggle"/treat/*.wig `;

}




