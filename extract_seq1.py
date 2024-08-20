#Goal: extract the DNA sequence of the ecDNA amplicon. Run extract_seq1.py followed by extract_seq2.sh.

#Part 0: Import packages.
import os
import pysam
import argparse
import pandas as pd

#Part 1: Starting from AA_output, grab all .bed files corresponding to the type of amplicon we want to analyze.
parser = argparse.ArgumentParser(description='Extract the DNA sequence of the amplicon.')
parser.add_argument('sample', metavar ='sample name', required = True, dest ='sample', action='store', 
                    help='Sample name. e.g. H143_rAsc-P4', type='string')
parser.add_argument('-d', '--dir', metavar ='directory', required = True, dest ='directory', action='store', 
                    help='Pass in the parent directory containing sample_output from AmpliconArchitect. e.g. home/H143_rAsc-P4', 
                    nargs=1, type='string')
parser.add_argument('--amp-type', required = False, dest ='amplicon-type', action='store', 
                    help='Type of amplicon to examine. e.g. BFB. ecDNA is the default.', 
                    nargs=1, default='ecDNA', type='string')
args = parser.parse_args()

sample = args.sample
parent_directory = args.dir
amp_type = args.amp-type
output_folder = sample + "_output/"
classification_folder = sample + "_classification/"

# get list of which amplicons match desired amp-type
path_to_result_table = os.path.join(parent_directory, output_folder, classification_folder, sample + "_result_table.tsv")
result_table_df = pd.read_csv(path_to_result_table, sep='\t')
result_table_df = result_table_df[result_table_df['Classification'] == amp_type]
filtered_amplicons = list(result_table_df.loc[:, 'AA amplicon number'])
filtered_amplicons.sort()

# acquire bed folder and bam files
path_to_bed_folder = os.path.join(parent_directory, output_folder, classification_folder, sample + "_classification_bed_files")
path_to_bam_file = os.path.join(parent_directory, output_folder, sample + ".cs.rmdup.bam")
imported = pysam.AlignmentFile(path_to_bam_file, mode = 'rb')

# create a new output folder
try:
    new_folder = os.path.join(parent_directory, sample + "_" + amp_type + "_sequences")
    os.mkdir(new_folder)
    print("Folder %s created!" % new_folder)
except FileExistsError:
    print("Folder %s already exists." % new_folder)
        
# extract the .bam reads from the .bed files of each amplicon
for i in filtered_amplicons:
    path_to_bed_file = os.path.join(path_to_bed_folder, sample + "_amplicon" + str(filtered_amplicons[i]) + "_" + amp_type + "_1_intervals.bed")
    curr_bed = pd.read_csv(path_to_bed_file, sep='\t', header=None)

    # create new folder and new files
    try:
        amplicon_folder = os.path.join(new_folder, sample + "_amplicon" + str(filtered_amplicons[i]))
        os.mkdir(amplicon_folder)
        print("Folder %s created!" % amplicon_folder)
    except FileExistsError:
        print("Folder %s already exists." % amplicon_folder)
    
    new_bam = os.path.join(amplicon_folder, sample + "_amplicon" + str(filtered_amplicons[i]) + ".bam")
    sorted_bam = os.path.join(amplicon_folder, sample + "_amplicon" + str(filtered_amplicons[i]) + "_sorted.bam")
    new_bam_index = sorted_bam + ".bai"
    new_fq1 = os.path.join(amplicon_folder, sample + "_amplicon" + str(filtered_amplicons[i]) + "_R1.fastq")
    new_fq2 = os.path.join(amplicon_folder, sample + "_amplicon" + str(filtered_amplicons[i]) + "_R2.fastq")

    # extract all reads corresponding to the regions identified in the bed file
    with open(new_bam, 'wb') as out_bam:
        bam_output = pysam.view("-b", "-h", "-M", "-L", path_to_bed_file, path_to_bam_file)
        out_bam.write(bam_output)

    pysam.sort("-o", sorted_bam, new_bam) #sort by genomic coordinates
    pysam.index(sorted_bam, new_bam_index)

    with open(new_fq1, 'wb') as fq1, open(new_fq2, 'wb') as fq2:
        pysam.fastq("-1", new_fq1, "-2", new_fq2, sorted_bam)


# my_bam_file = '/path/to/my/bam_file.bam'
# imported = pysam.AlignmentFile(my_bam_file, mode = 'rb')
# regions = ('1:2010000-20200000','2:2010000-20200000')
# alignments = []
# for region in regions:
#     bam = imported.fetch(region = region, until_eof = True)
#     alignments.extend([alignment for alignment in bam])




