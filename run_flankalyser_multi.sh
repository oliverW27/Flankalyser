#!/bin/bash

FLANKS_DIR=$1
FLANK_LEN=$2
OUTDIR=$3
for file in $FLANKS_DIR/*.txt; do
    [ -e "$file" ] || continue
    base=$(basename "$file")
    python3 flankalyser.py --mge_name "${base}" --summary "$file" --flank_len "$FLANK_LEN" --outdir "$OUTDIR"
done
