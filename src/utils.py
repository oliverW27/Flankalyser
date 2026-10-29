import logging
import random
import sys
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(seq):
    """Reverse complement of an uppercase ACGT sequence."""
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def parse_summary_file(summary_file):
    """Read the input TSV (acc_name, left_flank, right_flank) into three lists.

    Fields are kept raw, so pipe-separated multi-flank entries and 'none'
    are left for each analysis to handle. Exits if the file is missing or
    unreadable.
    """
    acc_names, l_data, r_data = [], [], []

    if not Path(summary_file).exists():
        logger.error(f"Summary file not found: {summary_file}")
        sys.exit(1)

    try:
        with open(summary_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    logger.warning(
                        f"Skipping malformed line {line_num} in {summary_file}: {line.strip()}"
                    )
                    continue
                acc_names.append(cols[0])
                l_data.append(cols[1])
                r_data.append(cols[2])
    except Exception as e:
        logger.exception(f"Failed while reading summary file: {e}")
        sys.exit(1)

    logger.info(
        f"Parsed {len(acc_names)} entries from {summary_file} "
        f"({len(l_data)} left, {len(r_data)} right flanks)"
    )
    return l_data, r_data, acc_names


def pair_flanks(left_data, right_data, flank_len):
    """Join left and right flanks from the same input line into insertion-site
    sequences (length 2 * flank_len), deduplicated, in order of first appearance.

    Pairing has to happen before deduplication, since the left/right
    correspondence on each line must be kept. Pipe-separated segments are
    paired by position. Segments of the wrong length or containing non-ACGT
    bases (including N and 'none') are skipped.
    """
    seen, skipped = {}, 0
    for left_field, right_field in zip(left_data, right_data):
        for left, right in zip(left_field.upper().split("|"), right_field.upper().split("|")):
            left, right = left.strip(), right.strip()
            combined = left + right
            if len(left) != flank_len or len(right) != flank_len or set(combined) - set("ACGT"):
                skipped += 1
                continue
            seen[combined] = None

    logger.info(
        f"Paired {len(seen)} unique insertion-site sequences of {2 * flank_len} bp "
        f"({skipped} segment(s) skipped: wrong length or non-ACGT)"
    )
    return list(seen)


def dinucleotide_shuffle(seq, max_attempts=200):
    """Shuffle a sequence while preserving its dinucleotide frequencies (Altschul & Erickson 1985).
    """
    seq = seq.upper()
    for _ in range(max_attempts):
        edges = defaultdict(list)
        for i in range(len(seq) - 1):
            edges[seq[i]].append(seq[i + 1])
        for base in edges:
            random.shuffle(edges[base])
        result = [seq[0]]
        current = seq[0]
        while edges[current]:
            next_base = edges[current].pop()
            result.append(next_base)
            current = next_base
        if len(result) == len(seq):
            return "".join(result)
    raise RuntimeError(f"Dinucleotide shuffle failed after {max_attempts} attempts")
