#!/bin/bash

split_reads() {
    local file="$1"

    if [[ ! -f "$file" ]]; then
        echo "File not found!"
        return 1
    fi

    local total_lines=$(wc -l < "$file")

    local half_lines=$((total_lines / 2))

    head -n "$half_lines" "$file" > ${1}_01.txt
    tail -n +$((half_lines + 1)) "$file" > ${1}_02.txt

    echo "File '${file}' has been split into '${file}_01.txt' and '${file}_02.txt"
}

split_reads "$1"