"""
Palindrome / inverted-repeat detection module.

Standalone script: reads a trimmed flanks TSV (acc_name\\tleft\\tright),
concatenates left+right per accession, searches for inverted repeats with
optional spacer and mismatches, compares against a dinucleotide-shuffled null
across a range of arm lengths, and writes:
  - a per-element summary line (inflection length, peak enrichment)
  - a full enrichment curve TSV
  - an enrichment ratio plot
  - a raw counts plot

Usage:
    python palindrome_analysis.py --input flanks.tsv --mge_name IS1397 --outdir results
"""

import argparse
import csv
import logging
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.utils import dinucleotide_shuffle

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


# ---------------------------------------------------------------------------
# Inverted repeat search
# ---------------------------------------------------------------------------

def _is_complement(a, b):
    return _COMPLEMENT.get(a) == b


def find_longest_inverted_repeat(seq, max_spacer, min_arm=3, max_mismatches=0,
                                 search_start=None, search_end=None,
                                 fixed_centre=None):
    """Find the longest inverted repeat (arms are reverse complements) with a spacer of 0+ nt.

    Optionally restrict the search to a window [search_start, search_end).
    If fixed_centre is set, only that index is tried as the palindrome centre.
    Returns dict with arm_len, left_start, right_end, spacer_size, arm_seq, mismatches.
    """
    seq = seq.upper()
    n = len(seq)
    if search_start is None:
        search_start = 0
    if search_end is None:
        search_end = n
    search_start = max(0, search_start)
    search_end = min(n, search_end)

    best = {"arm_len": 0, "left_start": None, "right_end": None,
            "spacer_size": None, "arm_seq": None, "mismatches": 0}

    centres = [fixed_centre] if fixed_centre is not None else range(search_start, search_end)
    for centre in centres:
        for spacer_size in range(0, max_spacer + 1):
            left_end = centre - spacer_size // 2
            right_start = centre + (spacer_size + 1) // 2
            if left_end < search_start or right_start >= search_end:
                break
            arm_len = 0
            mismatches = 0
            l, r = left_end - 1, right_start
            while l >= search_start and r < search_end:
                if _is_complement(seq[l], seq[r]):
                    arm_len += 1
                elif mismatches < max_mismatches:
                    mismatches += 1
                    arm_len += 1
                else:
                    break
                l -= 1
                r += 1
            if arm_len >= min_arm and arm_len > best["arm_len"]:
                best = {
                    "arm_len": arm_len,
                    "left_start": l + 1,
                    "right_end": r,
                    "spacer_size": spacer_size,
                    "arm_seq": seq[l + 1:r],
                    "mismatches": mismatches,
                }
    return best


# ---------------------------------------------------------------------------
# Input loading
# ---------------------------------------------------------------------------

def load_sequences(path, min_len=20):
    """Load flanks from TSV, pair L+R positionally, concatenate, deduplicate.

    Handles both single-flank rows (acc  left  right) and pipe-separated
    rows (acc  L1|L2|...  R1|R2|...), pairing by position.
    Returns list of unique concatenated L+R sequences.
    """
    seen = set()
    with open(path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue

            lefts = [s.strip() for s in row[1].split("|")]
            rights = [s.strip() for s in row[2].split("|")]

            for i in range(min(len(lefts), len(rights))):
                lf, rf = lefts[i], rights[i]

                # Skip empty, None, or all-N
                if not lf or not rf:
                    continue
                if lf.upper() == 'NONE' or rf.upper() == 'NONE':
                    continue

                lf_clean = lf.strip("Nn")
                rf_clean = rf.rstrip("Nn")

                if not lf_clean or not rf_clean:
                    continue
                if "N" in lf_clean.upper() or "N" in rf_clean.upper():
                    continue

                concat = (lf_clean + rf_clean).upper()
                if len(concat) < min_len:
                    continue

                seen.add(concat)

    return list(seen)


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def get_max_arm_lengths(seqs, max_spacer, max_mismatches, max_arm, centered=True):
    """For each sequence, find the longest palindrome at the insertion site.
    
    Returns list of arm lengths (0 if no palindrome >= 3bp found).
    """
    arm_lengths = []
    for s in seqs:
        fixed_centre = len(s) // 2 if centered else None
        res = find_longest_inverted_repeat(
            s,
            max_spacer=max_spacer,
            min_arm=3,
            max_mismatches=max_mismatches,
            fixed_centre=fixed_centre,
        )
        arm_lengths.append(min(res["arm_len"], max_arm))
    return arm_lengths


def compute_enrichment_curve(real_arms, shuffled_arms_list, min_arm, max_arm):
    """Compute enrichment ratio at each arm length threshold.
    
    At each threshold, counts how many sequences have a palindrome
    with arm length >= that threshold, for real vs shuffled.
    """
    curve = []
    for threshold in range(min_arm, max_arm + 1):
        real_count = sum(1 for a in real_arms if a >= threshold)

        shuf_counts = [sum(1 for a in shuf_arms if a >= threshold)
                       for shuf_arms in shuffled_arms_list]
        shuf_mean = np.mean(shuf_counts) if shuf_counts else 0.0
        shuf_sd = np.std(shuf_counts) if shuf_counts else 0.0

        if shuf_mean > 0:
            ratio = real_count / shuf_mean
        elif real_count > 0:
            ratio = float('inf')
        else:
            ratio = 1.0

        curve.append({
            'arm_length': threshold,
            'real_count': real_count,
            'shuffled_mean': shuf_mean,
            'shuffled_sd': shuf_sd,
            'ratio': ratio,
        })
    return curve


def extract_summary(curve, threshold=2.0):
    """Extract inflection length and peak enrichment from the curve."""
    inflection = None
    peak_enrichment = 0.0
    peak_at_length = None

    for point in curve:
        ratio = point['ratio']
        arm = point['arm_length']

        if ratio != float('inf') and ratio > peak_enrichment:
            peak_enrichment = ratio
            peak_at_length = arm

        if inflection is None and ratio >= threshold:
            inflection = arm

    # Handle inf cases
    if any(p['ratio'] == float('inf') for p in curve):
        inf_points = [p for p in curve if p['ratio'] == float('inf')]
        if not peak_at_length or peak_enrichment < 100:
            peak_enrichment = float('inf')
            peak_at_length = inf_points[0]['arm_length']

    return {
        'inflection_length': inflection,
        'peak_enrichment': peak_enrichment,
        'peak_at_length': peak_at_length,
    }


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _sig_stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_palindrome_analysis(curve, mge_name, summary, threshold, output_path):
    """Combined two-panel figure: raw counts (top) and enrichment ratio (bottom)."""
    arms = [p['arm_length'] for p in curve]
    ratios = [p['ratio'] for p in curve]
    real_counts = [p['real_count'] for p in curve]
    shuf_means = [p['shuffled_mean'] for p in curve]
    shuf_sds = [p['shuffled_sd'] for p in curve]

    # Cap inf for plotting
    max_finite = max((r for r in ratios if r != float('inf')), default=10)
    plot_cap = max_finite * 1.5
    plot_ratios = [min(r, plot_cap) if r != float('inf') else plot_cap
                   for r in ratios]

    x = np.arange(len(arms))
    width = 0.35

    fig, (ax_counts, ax_enrich) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # --- Top panel: raw counts ---
    ax_counts.bar(x - width/2, real_counts, width, label='Real', color='#2166ac')
    ax_counts.bar(x + width/2, shuf_means, width, label='Shuffled (mean)',
                  color='tab:orange', yerr=shuf_sds, capsize=2)
    ax_counts.set_ylabel('Sequences with palindrome ≥ arm length')
    ax_counts.set_title(f'{mge_name}: palindrome counts and enrichment by arm length')
    ax_counts.legend(loc='upper right', fontsize=8)
    ax_counts.grid(True, axis='y', alpha=0.3)

    # --- Bottom panel: enrichment ratio ---
    ax_enrich.bar(x, plot_ratios, width * 2, color='#2166ac', alpha=0.8, edgecolor='white')

    ax_enrich.axhline(1.0, color='grey', linestyle='--', linewidth=1, alpha=0.6,
                      label='No enrichment')
    ax_enrich.axhline(threshold, color='red', linestyle='--', linewidth=1, alpha=0.6,
                      label=f'Threshold ({threshold}x)')

    if summary['inflection_length'] is not None:
        infl_idx = arms.index(summary['inflection_length'])
        ax_enrich.axvline(infl_idx, color='green', linestyle='--',
                          linewidth=1.5, alpha=0.7,
                          label=f'Inflection (arm = {summary["inflection_length"]})')
        # Also mark on top panel
        ax_counts.axvline(infl_idx, color='green', linestyle='--',
                          linewidth=1.5, alpha=0.4)

    for i, r in enumerate(ratios):
        if r == float('inf'):
            ax_enrich.annotate('∞', xy=(i, plot_ratios[i]),
                               ha='center', va='bottom', fontsize=10, color='red')

    ax_enrich.set_xticks(x)
    ax_enrich.set_xticklabels(arms)
    ax_enrich.set_xlabel('Minimum arm length (bp)')
    ax_enrich.set_ylabel('Enrichment (real / shuffled)')
    ax_enrich.legend(loc='upper left', fontsize=8)
    ax_enrich.grid(True, axis='y', alpha=0.3)

    # Annotation
    infl_str = str(summary['inflection_length']) if summary['inflection_length'] else 'none'
    peak_str = (f"{summary['peak_enrichment']:.1f}x"
                if summary['peak_enrichment'] != float('inf') else '∞')
    peak_at = str(summary['peak_at_length']) if summary['peak_at_length'] else '—'
    ax_enrich.text(0.98, 0.95,
                   f"Inflection: {infl_str} bp\nPeak: {peak_str} at {peak_at} bp",
                   transform=ax_enrich.transAxes, ha='right', va='top', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.85))

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_palindrome_analysis(
    input_tsv, mge_name, output_dir,
    min_arm=3, max_arm=25, max_spacer=2, max_mismatches=0,
    n_shuffles=100, threshold=2.0, seed=42,
    min_seq_len=20, centered=True,
):
    """Run palindrome enrichment analysis across arm lengths."""
    logger = logging.getLogger("palindrome")

    if seed is not None:
        random.seed(seed)

    seqs = load_sequences(input_tsv, min_len=min_seq_len)
    if not seqs:
        logger.error(f"No valid sequences loaded from {input_tsv}")
        return None

    seq_len = len(seqs[0])
    logger.info(f"Loaded {len(seqs)} unique sequences ({seq_len} bp)")
    logger.info(
        f"Settings: arm {min_arm}-{max_arm}, spacer <= {max_spacer}, "
        f"mismatches <= {max_mismatches}, {n_shuffles} shuffles, "
        f"centred={'yes' if centered else 'no'}"
    )

    # --- Real sequences ---
    logger.info("Scanning real sequences...")
    real_arms = get_max_arm_lengths(seqs, max_spacer, max_mismatches, max_arm, centered)

    # --- Shuffled controls ---
    logger.info(f"Scanning {n_shuffles} shuffled replicates...")
    shuffled_arms_list = []
    for i in range(n_shuffles):
        shuffled_seqs = []
        for s in seqs:
            try:
                shuffled_seqs.append(dinucleotide_shuffle(s))
            except RuntimeError:
                shuffled_seqs.append(s)
        shuf_arms = get_max_arm_lengths(shuffled_seqs, max_spacer, max_mismatches,
                                         max_arm, centered)
        shuffled_arms_list.append(shuf_arms)
        if (i + 1) % 10 == 0:
            logger.info(f"  {i+1}/{n_shuffles} shuffles done")

    # --- Enrichment curve ---
    curve = compute_enrichment_curve(real_arms, shuffled_arms_list, min_arm, max_arm)
    summary = extract_summary(curve, threshold)

    summary['n_seqs'] = len(seqs)
    summary['pct_at_inflection'] = None
    if summary['inflection_length'] is not None:
        real_count_at_infl = next((p['real_count'] for p in curve
                                    if p['arm_length'] == summary['inflection_length']), None)
        if real_count_at_infl is not None:
            summary['pct_at_inflection'] = real_count_at_infl / len(seqs) * 100

    peak_str = (f"{summary['peak_enrichment']:.1f}x"
                if summary['peak_enrichment'] != float('inf') else 'inf')
    logger.info(
        f"Result: inflection={summary['inflection_length']}, "
        f"peak={peak_str} at {summary['peak_at_length']}bp"
    )

    # --- Output ---
    out_dir = Path(output_dir)
    files_dir = out_dir / "files"
    graphs_dir = out_dir / "graphs"
    files_dir.mkdir(parents=True, exist_ok=True)
    graphs_dir.mkdir(parents=True, exist_ok=True)

    # Parameter suffix so different parameter combos never overwrite each other
    param_suffix = f"arm{min_arm}_mm{max_mismatches}_sp{max_spacer}"

    # Summary line (appends to running file, one file per parameter combo)
    summary_path = files_dir / f"palindrome_enrichment_summary_{param_suffix}.tsv"
    write_header = not summary_path.exists()
    with open(summary_path, "a") as f:
        if write_header:
            f.write("element\tn_seqs\tinflection_length\tpeak_enrichment\tpeak_at_length\n")
        infl_str = str(summary['inflection_length']) if summary['inflection_length'] else 'none'
        f.write(f"{mge_name}\t{len(seqs)}\t{infl_str}\t{peak_str}\t"
                f"{summary['peak_at_length']}\n")
    logger.info(f"Summary appended to {summary_path}")

    # Full curve
    curve_path = files_dir / f"{mge_name}_PALINDROME_ENRICHMENT_CURVE_{param_suffix}.tsv"
    with open(curve_path, "w") as f:
        f.write("element\tarm_length\treal_count\tshuffled_mean\tshuffled_sd\tratio\n")
        for p in curve:
            ratio_str = f"{p['ratio']:.2f}" if p['ratio'] != float('inf') else 'inf'
            f.write(f"{mge_name}\t{p['arm_length']}\t{p['real_count']}\t"
                    f"{p['shuffled_mean']:.2f}\t{p['shuffled_sd']:.2f}\t{ratio_str}\n")
    logger.info(f"Enrichment curve written to {curve_path}")

    # Plot (combined counts + enrichment)
    plot_path = graphs_dir / f"{mge_name}_PALINDROME_ENRICHMENT_{param_suffix}.png"
    plot_palindrome_analysis(curve, mge_name, summary, threshold, plot_path)
    logger.info(f"Palindrome plot written to {plot_path}")

    logger.info("Palindrome analysis complete.")
    return summary


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="Palindromic enrichment by arm length, compared to dinucleotide-shuffled null."
    )
    p.add_argument("--input", "-i", required=True,
                   help="Trimmed flanks TSV: acc_name<TAB>left<TAB>right")
    p.add_argument("--mge_name", "-m", required=True,
                   help="MGE name used in output filenames")
    p.add_argument("--outdir", "-o", default=".",
                   help="Output directory (default: current)")
    p.add_argument("--min_arm", type=int, default=3,
                   help="Minimum arm length to report (default: 3)")
    p.add_argument("--max_arm", type=int, default=25,
                   help="Maximum arm length to scan (default: 25)")
    p.add_argument("--max_spacer", type=int, default=2,
                   help="Maximum spacer size (default: 2)")
    p.add_argument("--max_mismatches", type=int, default=0,
                   help="Maximum mismatches per arm (default: 0)")
    p.add_argument("--n_shuffles", type=int, default=100,
                   help="Number of null shuffles (default: 100)")
    p.add_argument("--threshold", type=float, default=2.0,
                   help="Enrichment ratio threshold for inflection (default: 2.0)")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed (default: 42)")
    p.add_argument("--min_seq_len", type=int, default=20,
                   help="Skip sequences shorter than this (default: 20)")
    p.add_argument("--no_centered", action="store_true",
                   help="Search full sequence instead of only insertion midpoint")
    return p.parse_args()


def main():
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stderr,
    )
    run_palindrome_analysis(
        input_tsv=args.input,
        mge_name=args.mge_name,
        output_dir=args.outdir,
        min_arm=args.min_arm,
        max_arm=args.max_arm,
        max_spacer=args.max_spacer,
        max_mismatches=args.max_mismatches,
        n_shuffles=args.n_shuffles,
        threshold=args.threshold,
        seed=args.seed,
        min_seq_len=args.min_seq_len,
        centered=not args.no_centered,
    )


if __name__ == "__main__":
    main()