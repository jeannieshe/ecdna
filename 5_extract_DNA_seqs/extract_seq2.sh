#!/bin/bash

# Check if the csv file argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

csv_file="$1"

# Check if the file exists
if [ ! -f "$csv_file" ]; then
    echo "File not found: $csv_file"
    exit 1
fi

#acquire .vcf from .fastq using freebayes
#ref.fa = 2024 genome.fa. --min-coverage 5 
freebayes -f ref2.fa --min-coverage 5 D327_r1PlfL1_amplicon2.bam >var.vcf


#look for the soft clipped reads from the bam file/vcf file (they align up to a certain point)
bcftools view -Ou --types=snps .vcf.gz | bcftools stats - | awk '/^SN/'

#check in cigar string which are soft clipped and which are not


#run unicycler to assemble the reads into a circular template
unicycler -1 D327_r1PlfL1_amplicon2_R1.fastq -2 D327_r1PlfL1_amplicon2_R2.fastq -o /home/jeannie/group/OvCa_Genomic_Data/WGS_merged/D327_r1PlfL1/D327_r1PlfL1_ecDNA_sequences/D327_r1PlfL1_amplicon2
#output: fasta
# try installing unicycler on the macos to see if i can run it on the computer.
#then try qemu for downloading a linux subsystem to run unicycler instead on linux.

#sort and align with minimap2.
# input: fasta, reference fasta
# output: sam (verified the fastq)