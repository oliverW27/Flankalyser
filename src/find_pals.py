"""
Palindrome / inverted-repeat enrichment analysis.

Pairs left+right flanks into insertion-site sequences, finds the longest
inverted repeat in each, and compares how many sequences have one at least
N bp long against a dinucleotide-shuffled null, for every N from MIN_ARM to
MAX_ARM. Writes the full enrichment curve and a plot of real vs. shuffled counts
by arm length.
"""

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from src.utils import COMPLEMENT, dinucleotide_shuffle, pair_flanks

logger = logging.getLogger(__name__)

MIN_ARM = 3    # shortest arm counted as a palindrome, and the start of the enrichment curve
MAX_ARM = 50   # longer arms are counted as this length
INF = float('inf')


def longest_arm(seq, max_spacer, max_mismatches):
    """Arm length of the longest inverted repeat in seq (arms are reverse
    complements, separated by a spacer of 0..max_spacer nt), trying every
    position as the centre. Returns 0 if none is at least MIN_ARM long."""
    n = len(seq)
    best = 0
    for centre in range(n):
        for spacer_size in range(max_spacer + 1):
            l = centre - spacer_size // 2 - 1
            r = centre + (spacer_size + 1) // 2
            if l < -1 or r >= n:
                break
            arm_len = mismatches = 0
            while l >= 0 and r < n:
                if COMPLEMENT.get(seq[l]) != seq[r]:
                    if mismatches == max_mismatches:
                        break
                    mismatches += 1
                arm_len += 1
                l -= 1
                r += 1
            best = max(best, arm_len)
    return best if best >= MIN_ARM else 0


def get_max_arm_lengths(seqs, max_spacer, max_mismatches):
    """Longest palindrome arm in each sequence, capped at MAX_ARM."""
    return [min(longest_arm(s, max_spacer, max_mismatches), MAX_ARM) for s in seqs]


def compute_enrichment_curve(real_arms, shuffled_arms_list):
    """For each arm length from MIN_ARM to MAX_ARM, count how many sequences
    have a palindrome at least that long, real vs. shuffled."""
    curve = []
    for arm_length in range(MIN_ARM, MAX_ARM + 1):
        real_count = sum(1 for a in real_arms if a >= arm_length)
        shuf_counts = [sum(1 for a in shuf_arms if a >= arm_length)
                       for shuf_arms in shuffled_arms_list]
        shuf_mean = np.mean(shuf_counts) if shuf_counts else 0.0
        shuf_sd = np.std(shuf_counts) if shuf_counts else 0.0

        # inf = real palindromes this long, but none in any shuffle
        if shuf_mean > 0:
            ratio = real_count / shuf_mean
        elif real_count > 0:
            ratio = INF
        else:
            ratio = 1.0

        curve.append({
            'arm_length': arm_length,
            'real_count': real_count,
            'shuffled_mean': shuf_mean,
            'shuffled_sd': shuf_sd,
            'ratio': ratio,
        })
    return curve


def extract_summary(curve, threshold):
    """Inflection = shortest arm length whose enrichment reaches threshold.
    Peak = highest enrichment and the arm length it occurs at."""
    inflection = next((p['arm_length'] for p in curve if p['ratio'] >= threshold), None)

    finite = [p for p in curve if p['ratio'] != INF]
    peak = max(finite, key=lambda p: p['ratio'], default=None)
    peak_enrichment = peak['ratio'] if peak and peak['ratio'] > 0 else 0.0
    peak_at_length = peak['arm_length'] if peak and peak['ratio'] > 0 else None

    # An infinite ratio is reported as the peak, unless a finite peak is already very large
    inf_arms = [p['arm_length'] for p in curve if p['ratio'] == INF]
    if inf_arms and peak_enrichment < 100:
        peak_enrichment, peak_at_length = INF, inf_arms[0]

    return {
        'inflection_length': inflection,
        'peak_enrichment': peak_enrichment,
        'peak_at_length': peak_at_length,
    }


def _fmt_ratio(r):
    return 'inf' if r == INF else f"{r:.1f}x"


def plot_palindrome_analysis(curve, mge_name, summary, output_path):
    """Bar chart of how many sequences have a palindrome at least each arm length,
    real vs. shuffled, with the inflection arm length marked."""
    arms = [p['arm_length'] for p in curve]
    x = np.arange(len(arms))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.bar(x - width/2, [p['real_count'] for p in curve], width, label='Real', color='#2166ac')
    ax.bar(x + width/2, [p['shuffled_mean'] for p in curve], width, label='Shuffled (mean)',
           color='tab:orange', yerr=[p['shuffled_sd'] for p in curve], capsize=2)

    infl = summary['inflection_length']
    if infl is not None:
        ax.axvline(arms.index(infl), color='green', linestyle='--', linewidth=1.5, alpha=0.4,
                   label=f'Inflection ({infl} bp)')

    ticks = [i for i, a in enumerate(arms) if i == 0 or a % 5 == 0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([arms[i] for i in ticks])
    ax.set_xlabel('Minimum arm length (bp)')
    ax.set_ylabel('Sequences with palindrome ≥ arm length')
    ax.set_title(f'{mge_name}: palindrome counts by arm length')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    plt.close()


def run_palindrome_analysis(left_data, right_data, flank_len, mge_name, output_dir,
                            max_spacer, max_mismatches, n_shuffles, threshold):
    """Palindrome enrichment pipeline: pair flanks, scan real and shuffled
    sequences, then write the curve TSV and plot.

    Args:
        left_data, right_data: raw per-line fields from the summary file
        threshold: enrichment ratio that defines the inflection arm length
    """
    seqs = pair_flanks(left_data, right_data, flank_len)
    if not seqs:
        logger.error("No usable left/right flank pairs -- skipping palindrome analysis")
        return

    logger.info(f"Loaded {len(seqs)} unique sequences ({len(seqs[0])} bp)")
    logger.info(
        f"Settings: arm {MIN_ARM}-{MAX_ARM}, spacer <= {max_spacer}, "
        f"mismatches <= {max_mismatches}, {n_shuffles} shuffles"
    )

    logger.info("Scanning real sequences...")
    real_arms = get_max_arm_lengths(seqs, max_spacer, max_mismatches)

    # Fall back to the real sequence if a shuffle fails (rare; see utils.dinucleotide_shuffle)
    logger.info(f"Scanning {n_shuffles} shuffled replicates...")
    shuffled_arms_list = []
    for i in range(n_shuffles):
        shuffled_seqs = []
        for s in seqs:
            try:
                shuffled_seqs.append(dinucleotide_shuffle(s))
            except RuntimeError:
                shuffled_seqs.append(s)
        shuffled_arms_list.append(get_max_arm_lengths(shuffled_seqs, max_spacer, max_mismatches))
        if (i + 1) % 10 == 0:
            logger.info(f"  {i+1}/{n_shuffles} shuffles done")

    curve = compute_enrichment_curve(real_arms, shuffled_arms_list)
    summary = extract_summary(curve, threshold)

    # % of sequences with a palindrome at least as long as the inflection arm length
    infl = summary['inflection_length']
    pct_str = f"{curve[infl - MIN_ARM]['real_count'] / len(seqs):.1%}" if infl is not None else "NA"
    logger.info(
        f"Result: inflection={infl}, % at inflection={pct_str}, "
        f"peak={_fmt_ratio(summary['peak_enrichment'])} at {summary['peak_at_length']}bp"
    )

    files_dir = Path(output_dir) / "files"
    graphs_dir = Path(output_dir) / "graphs"
    # Different parameter combinations get their own files, so they never overwrite each other
    param_suffix = f"mm{max_mismatches}_sp{max_spacer}"

    curve_path = files_dir / f"{mge_name}_PALINDROME_ENRICHMENT_CURVE_{param_suffix}.tsv"
    with open(curve_path, "w") as f:
        f.write("arm_length\treal_count\tshuffled_mean\tshuffled_sd\tratio\n")
        for p in curve:
            if p['real_count'] == 0 and p['shuffled_mean'] == 0:
                ratio_str = 'NA'  # no palindromes this long in real or shuffled sequences
            elif p['ratio'] == INF:
                ratio_str = 'inf'
            else:
                ratio_str = f"{p['ratio']:.2f}"
            f.write(f"{p['arm_length']}\t{p['real_count']}\t"
                    f"{p['shuffled_mean']:.2f}\t{p['shuffled_sd']:.2f}\t{ratio_str}\n")
    logger.info(f"Enrichment curve written to {curve_path}")

    plot_path = graphs_dir / f"{mge_name}_PALINDROME_COUNTS_{param_suffix}.png"
    plot_palindrome_analysis(curve, mge_name, summary, plot_path)
    logger.info(f"Palindrome plot written to {plot_path}")
