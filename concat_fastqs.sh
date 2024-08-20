#!/bin/bash

concat_fq_paired() {
    if [[ "$#" -lt 2 ]]; then
        echo "Incorrect number of arguments passed. Usage: ./concat_fastqs.sh <directory> <sample name> [-o <output directory if different from input>]"
        return 1
    fi

    if [[ ! -d "$1" ]]; then
        echo "Directory not found"
        return 1
    fi

    dir="$1"
    sample="$2"
    dest="$dir" #set default destination directory

    echo "Input directory $1."

    shift 2
    while getopts "o:" flag; do
        echo "Reviewing the flag."
        case $flag in
            o) #if an output directory is passed
                output="$OPTARG"
                if [[ "$output" == "" ]]; then
                    dest="$dir"
                    echo "No new directory made. Output directory set to input directory."
                fi
                if [ ! -d "$output"  ]; then
                    mkdir -p "$output"
                    echo "New directory $output made, set to output directory."
                    dest="$output"
                fi
                if [ -d "$output"  ]; then
                    echo "Output directory set."
                    dest="$output"
                fi
                ;;
            *) #unknown flag
                echo "Invalid option: -$OPTARG" >&2 #this specifically sends the message to stderr
                return 1
                ;;
        esac
    done

    files1=($(find "$dir" -type f -name '*_1.fq.gz' | sort ) )
    files2=()

    for i in "${files1[@]}"; do #remove the 1, replace with 2
        new="${i//_1.fq.gz/_2.fq.gz}" #substring replacement
        files2+=("$new")
    done

    # for f in "${files1[@]}"; do #if the last line is not a new line character
    #     ( cat "${f}"; echo ) >> "${sample}_1_merged.fq.gz";
    # done

    # for f in "${files2[@]}"; do
    #     ( cat "${f}"; echo ) >> "${sample}_2_merged.fq.gz";
    # done

    echo "Output directory $dest."

    cat "${files1[@]}" > "${dest}/${sample}_R1_merged.fq.gz"
    cat "${files2[@]}" > "${dest}/${sample}_R2_merged.fq.gz"
    echo "Files have been merged into '${dest}/${sample}_R1_merged.fq.gz' and '${dest}/${sample}_R2_merged.fq.gz'."
}

concat_fq_paired "$@"