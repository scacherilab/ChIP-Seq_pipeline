Dependencies:
1) The following Perl packages:

List::Util
Getopt::Usage
Cwd
IPC::System::Simple

2) fastx_toolkit (tested on v0.0.13)
3) bowtie, bowtie2 or bwa (depending on the aligner you use)
4) samtools (tested on v0.1.18)
5) macs2
6) wigToBigWig (for bigWIG generation from WIGs <--not necessary if --nowig option is selected)

For paired end samples, the 2 paired end files must be combined as processed as single end experiments. Paired-end option will be added in the next version
