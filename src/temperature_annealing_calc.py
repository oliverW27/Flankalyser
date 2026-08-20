#!/usr/bin/env python3
"""
Melting analysis of IS element flanking sequences.

Two per-position profiles centered on the insertion site, each shown against
a dinucleotide-shuffled control:
  1. Mean hydrogen bond count (2 = AT, 3 = GC)
  2. Nearest-neighbor ΔG (SantaLucia 1998 unified params)

Usage: python melting_analysis.py input.tsv [output_prefix]

Input: tab-separated (sample_id  left_flank  right_flank)
  Oriented as: LEFT_FLANK | IS | RIGHT_FLANK
  Ns are skipped; all-N sequences excluded; flanks deduplicated.

References:
  SantaLucia J Jr (1998) PNAS 95:1460-1465  (NN params)
"""

import logging
import random
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from src.utils import dinucleotide_shuffle

logger = logging.getLogger(__name__)

# --- SantaLucia (1998) unified NN ΔG°37 in kcal/mol ---
NN_DG = {
    'AA': -1.00, 'TT': -1.00, 'AT': -0.88, 'TA': -0.58,
    'CA': -1.45, 'TG': -1.45, 'GT': -1.44, 'AC': -1.44,
    'CT': -1.28, 'AG': -1.28, 'GA': -1.30, 'TC': -1.30,
    'CG': -2.17, 'GC': -2.24, 'GG': -1.84, 'CC': -1.84,
}
HBONDS = {'A': 2, 'T': 2, 'G': 3, 'C': 3}


def parse_input(filepath):
    """Read TSV, split | flanks, deduplicate, skip all-N."""
    lefts, rights = set(), set()
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            # Split pipe-separated flanks into independent sequences
            for lf in parts[1].upper().split('|'):
                lf = lf.strip()
                if lf and set(lf) != {'N'}:
                    lefts.add(lf)
            for rf in parts[2].upper().split('|'):
                rf = rf.strip()
                if rf and set(rf) != {'N'}:
                    rights.add(rf)
    return list(lefts), list(rights)


def profile_with_sem(seqs, per_seq_fn):
    """Mean and SEM across per-sequence profiles."""
    profiles = [per_seq_fn(s) for s in seqs if per_seq_fn(s) is not None]
    if not profiles:
        return None, None
    max_len = max(len(p) for p in profiles)
    arr = np.full((len(profiles), max_len), np.nan)
    for i, p in enumerate(profiles):
        arr[i, :len(p)] = p
    mean = np.nanmean(arr, axis=0)
    sem = np.nanstd(arr, axis=0) / np.sqrt(np.sum(~np.isnan(arr), axis=0))
    return mean, sem


def hbond_per_seq(seq):
    return np.array([HBONDS.get(b, np.nan) for b in seq])


def dg_per_seq(seq):
    vals = [NN_DG.get(seq[i:i+2], np.nan) for i in range(len(seq) - 1)]
    return np.array(vals)


def run_melting_analysis(left_seqs, right_seqs, mge_name, output_dir, seed=42):
    """Compute and plot melting profiles for a set of flank sequences.

    Two per-position profiles centered on the insertion site, each shown
    against a dinucleotide-shuffled control:
      1. Mean hydrogen bond count (2 = AT, 3 = GC)
      2. Nearest-neighbor ΔG (SantaLucia 1998 unified params)

    Args:
        left_seqs, right_seqs: lists of equal-length flank sequences (already
            deduplicated/filtered -- e.g. flankalyser.py's extract_unique_sequences)
        mge_name: name of the MGE, used as the output filename prefix
        output_dir: directory to save the plot under (writes to <output_dir>/graphs/)
        seed: random seed for the dinucleotide-shuffled control (for reproducibility)

    Returns:
        dict with mean H-bond count per side, or None if no usable sequences were found.
    """
    lefts = [s.upper() for s in left_seqs if 'N' not in s.upper()]
    rights = [s.upper() for s in right_seqs if 'N' not in s.upper()]
    logger.info(
        f"Melting analysis input: {len(lefts)}/{len(left_seqs)} left and "
        f"{len(rights)}/{len(right_seqs)} right flanks usable (N-free)"
    )
    if not lefts or not rights:
        logger.warning("No N-free sequences on one or both sides -- skipping melting analysis")
        return None

    if seed is not None:
        random.seed(seed)
    shuf_lefts = [dinucleotide_shuffle(s) for s in lefts]
    shuf_rights = [dinucleotide_shuffle(s) for s in rights]

    left_len = len(lefts[0])
    right_len = len(rights[0])

    # --- Compute all profiles ---
    logger.debug("Computing H-bond profiles...")
    hb_l, hb_l_sem = profile_with_sem(lefts, hbond_per_seq)
    hb_r, hb_r_sem = profile_with_sem(rights, hbond_per_seq)
    hb_sl, _ = profile_with_sem(shuf_lefts, hbond_per_seq)
    hb_sr, _ = profile_with_sem(shuf_rights, hbond_per_seq)

    logger.debug("Computing nearest-neighbor ΔG profiles...")
    dg_l, dg_l_sem = profile_with_sem(lefts, dg_per_seq)
    dg_r, dg_r_sem = profile_with_sem(rights, dg_per_seq)
    dg_sl, _ = profile_with_sem(shuf_lefts, dg_per_seq)
    dg_sr, _ = profile_with_sem(shuf_rights, dg_per_seq)

    # --- X-axes relative to insertion ---
    x_l = np.arange(-left_len, 0)
    x_r = np.arange(1, right_len + 1)
    x_l_dg = x_l[:-1] + 0.5   # dinucleotide step centers
    x_r_dg = x_r[:-1] + 0.5

    # --- Plot ---
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    dk = dict(linewidth=1.0)
    sk = dict(linewidth=1.0, linestyle='--', alpha=0.4, color='grey')
    fa = 0.12

    # Panel 1: H-bonds
    ax = axes[0]
    ax.plot(x_l, hb_l, color='#2166ac', **dk, label='Left flank')
    ax.fill_between(x_l, hb_l - hb_l_sem, hb_l + hb_l_sem, color='#2166ac', alpha=fa)
    ax.plot(x_r, hb_r, color='#b2182b', **dk, label='Right flank')
    ax.fill_between(x_r, hb_r - hb_r_sem, hb_r + hb_r_sem, color='#b2182b', alpha=fa)
    if hb_sl is not None:
        ax.plot(x_l, hb_sl, **sk, label='Shuffled')
        ax.plot(x_r, hb_sr, **sk)
    ax.axvline(0.5, color='grey', ls='--', alpha=0.4)
    ax.set_ylabel('Mean H-bonds per bp')
    ax.set_title('Hydrogen bond count (2 = AT, 3 = GC)')
    ax.legend(fontsize=8, loc='upper right')

    # Panel 2: NN ΔG  (less negative = easier to melt = higher on plot)
    ax = axes[1]
    ax.plot(x_l_dg, dg_l, color='#2166ac', **dk)
    ax.fill_between(x_l_dg, dg_l - dg_l_sem, dg_l + dg_l_sem, color='#2166ac', alpha=fa)
    ax.plot(x_r_dg, dg_r, color='#b2182b', **dk)
    ax.fill_between(x_r_dg, dg_r - dg_r_sem, dg_r + dg_r_sem, color='#b2182b', alpha=fa)
    if dg_sl is not None:
        ax.plot(x_l_dg, dg_sl, **sk)
        ax.plot(x_r_dg, dg_sr, **sk)
    ax.axvline(0.5, color='grey', ls='--', alpha=0.4)
    ax.set_ylabel('ΔG°₃₇ (kcal/mol)')
    ax.set_title('Nearest-neighbor stacking ΔG (SantaLucia 1998)')
    ax.annotate('↑ easier to melt', xy=(0.01, 0.95), xycoords='axes fraction',
                fontsize=7, color='grey', va='top')
    ax.set_xlabel(f'Position relative to {mge_name} insertion site')

    plt.tight_layout()
    graphs_dir = Path(output_dir) / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)
    png_path = graphs_dir / f"{mge_name}_MELTING_PROFILE.png"
    plt.savefig(png_path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Melting profile plot written to {png_path}")

    result = {
        "mean_hbonds_left": float(np.nanmean(hb_l)) if hb_l is not None else float("nan"),
        "mean_hbonds_right": float(np.nanmean(hb_r)) if hb_r is not None else float("nan"),
    }
    logger.info(
        f"Melting result: mean H-bonds left/right = {result['mean_hbonds_left']:.2f}/"
        f"{result['mean_hbonds_right']:.2f}"
    )
    return result


def main():
    """Standalone CLI: python temperature_annealing_calc.py input.tsv [output_prefix]."""
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s', stream=sys.stderr)

    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "melting_analysis"

    logger.info(f"Parsing {infile}...")
    lefts, rights = parse_input(infile)
    logger.info(f"{len(lefts)} unique left, {len(rights)} unique right flanks")

    run_melting_analysis(lefts, rights, prefix, ".")


if __name__ == '__main__':
    main()
