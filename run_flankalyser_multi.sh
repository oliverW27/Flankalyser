#!/bin/bash

FLANKS_DIR=$1
FLANK_LEN=$2

for file in $FLANKS_DIR/flanks_ISFinder_*.txt; do
    # echo "$FLANK_DIR"
    MGE=$(echo "$file" | sed -E 's/.*_(IS[0-9A-Za-z\-]+)(.txt)$/\1/')
    #echo "$MGE"
    python3 flankalyser.py --mge_name "${MGE}" --summary "$file" --flank_len "$FLANK_LEN" --outdir "$FLANKS_DIR" --adjust_proportions --gc_file expected_gc_per_element.tsv

done
