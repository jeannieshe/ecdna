#!/bin/bash
#in order to run this function, you should source the script: source clone_calling.sh. then you can call clone_calling with directories
clone_calling() {
    #this function creates a new CSV file and calls bam_parsing.R
    #the directories passed in as arguments should have the sample name as the folder name, and either "RNA" or "ATAC" in the name, e.g. /home/../H180_1_snATACseq
    #Usage: ./clone_calling.sh <dir> <dir> <dir>

    > bam_config.csv #clears the file

    echo "sample_id,bam_file,cell_barcode,assay,celltag_version" >> bam_config.csv #adds top line

    for i in "$@"; do
        dir="$i"
        if [[ ! -d $dir ]]; then
        echo "Directory not found"
        return 1
        fi
        
        sample="${echo $dir | rev | cut -d/ -f1 | rev}"

        if [[ $sample == *"RNA"* ]]; then
            echo "${sample},${dir}/outs/possorted_genome_bam.bam,${dir}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz,RNA,custom"; else
            echo "${sample},${dir}/outs/possorted_bam.bam,${dir}/outs/filtered_peak_bc_matrix/barcodes.tsv,ATAC,custom"
        fi
    done

}

clone_calling "$@"
cloneCalling_scripts/bam_parsing.R "$(pwd)/bam_config.csv"

# then call python3 celltag_analysis_single_assay_py.py after editing the appropriate sample names and sample type (lines 25 and 79)
