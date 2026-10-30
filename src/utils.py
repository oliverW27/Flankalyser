import logging
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
VALID_FLANK_CHARS = set("ACGTN|")  # '|' separates multiple flanks in one field


def reverse_complement(seq):
    """Reverse complement of an uppercase ACGT sequence."""
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def _clean_flank_field(field, tally):
    """Uppercase a flank field and turn anything other than A/C/G/T/N (e.g. ambiguity
    codes or gaps) into N, so every analysis only ever sees A/C/G/T/N."""
    tally["lowercase"] += sum(c.islower() for c in field)
    field = field.upper()
    odd = [c for c in field if c not in VALID_FLANK_CHARS]
    if odd:
        tally["converted"] += len(odd)
        tally["chars"].update(odd)
        field = "".join(c if c in VALID_FLANK_CHARS else "N" for c in field)
    return field


def parse_summary_file(summary_file):
    """Read the input TSV (acc_name, left_flank, right_flank) into three lists.

    Flanks are uppercased and any character other than A/C/G/T/N is replaced
    with N. Pipe-separated multi-flank entries are left for each analysis to
    handle. Exits if the file is missing or unreadable.
    """
    acc_names, l_data, r_data = [], [], []
    tally = {"lowercase": 0, "converted": 0, "chars": Counter()}

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
                l_data.append(_clean_flank_field(cols[1], tally))
                r_data.append(_clean_flank_field(cols[2], tally))
    except Exception as e:
        logger.exception(f"Failed while reading summary file: {e}")
        sys.exit(1)

    logger.info(
        f"Parsed {len(acc_names)} entries from {summary_file} "
        f"({len(l_data)} left, {len(r_data)} right flanks)"
    )
    if tally["lowercase"]:
        logger.info(f"Converted {tally['lowercase']:,} lowercase bases to uppercase")
    if tally["converted"]:
        chars = ", ".join(f"'{c}' x{n:,}" for c, n in tally["chars"].most_common())
        logger.warning(
            f"Replaced {tally['converted']:,} characters other than A/C/G/T/N with N ({chars})"
        )
    return l_data, r_data, acc_names


def pair_flanks(left_data, right_data, flank_len):
    """Join left and right flanks from the same input line into insertion-site
    sequences (length 2 * flank_len), deduplicated, in order of first appearance.

    Pairing has to happen before deduplication, since the left/right
    correspondence on each line must be kept. Pipe-separated segments are
    paired by position. Segments of the wrong length or containing non-ACGT
    bases (including N) are skipped.
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
