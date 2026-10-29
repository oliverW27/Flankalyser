"""
Melting analysis of flanking sequences.

Per-position nearest-neighbour stacking ΔG (SantaLucia 1998 unified
parameters) either side of the insertion site, shown against a
dinucleotide-shuffled control. Less negative ΔG = easier to melt.

References:
  SantaLucia J Jr (1998) PNAS 95:1460-1465
"""

import logging
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.utils import dinucleotide_shuffle

logger = logging.getLogger(__name__)

# SantaLucia (1998) unified nearest-neighbour ΔG°37, kcal/mol
NN_DG = {
    'AA': -1.00, 'TT': -1.00, 'AT': -0.88, 'TA': -0.58,
    'CA': -1.45, 'TG': -1.45, 'GT': -1.44, 'AC': -1.44,
    'CT': -1.28, 'AG': -1.28, 'GA': -1.30, 'TC': -1.30,
    'CG': -2.17, 'GC': -2.24, 'GG': -1.84, 'CC': -1.84,
}


def dg_per_seq(seq):
    """ΔG of each dinucleotide step in seq (NaN for steps with non-ACGT bases)."""
    return np.array([NN_DG.get(seq[i:i+2], np.nan) for i in range(len(seq) - 1)])


def profile_with_sem(seqs, per_seq_fn):
    """Per-position mean and SEM across equal-length sequences, ignoring NaNs."""
    arr = np.array([per_seq_fn(s) for s in seqs])
    n = np.sum(~np.isnan(arr), axis=0)
    return np.nanmean(arr, axis=0), np.nanstd(arr, axis=0) / np.sqrt(n)

def _shuffle_all(seqs):
    """Dinucleotide-shuffle each sequence, dropping any that fail to shuffle."""
    shuffled = []
    for s in seqs:
        try:
            shuffled.append(dinucleotide_shuffle(s))
        except RuntimeError:
            pass
    return shuffled


def run_melting_analysis(left_seqs, right_seqs, mge_name, output_dir):
    """Plot the mean ΔG profile of each flank against a dinucleotide-shuffled
    control, to <output_dir>/graphs/<mge_name>_MELTING_PROFILE.png.

    Args:
        left_seqs, right_seqs: equal-length unique flank sequences, from
            flankalyser.py's extract_unique_sequences
        mge_name: name of the MGE, used in the title and output filename
    """
    # Flanks containing any N are left out entirely
    lefts = [s.upper() for s in left_seqs if 'N' not in s.upper()]
    rights = [s.upper() for s in right_seqs if 'N' not in s.upper()]
    logger.info(
        f"Melting analysis input: {len(lefts)}/{len(left_seqs)} left and "
        f"{len(rights)}/{len(right_seqs)} right flanks usable (N-free)"
    )
    if not lefts or not rights:
        logger.warning("No N-free sequences on one or both sides -- skipping melting analysis")
        return

    dg_l, dg_l_sem = profile_with_sem(lefts, dg_per_seq)
    dg_r, dg_r_sem = profile_with_sem(rights, dg_per_seq)
    dg_sl, _ = profile_with_sem(_shuffle_all(lefts), dg_per_seq)
    dg_sr, _ = profile_with_sem(_shuffle_all(rights), dg_per_seq)

    # Bases are at -flank_len..-1 and +1..+flank_len, so each step sits halfway between two bases
    x_l = np.arange(-len(lefts[0]), -1) + 0.5
    x_r = np.arange(1, len(rights[0])) + 0.5

    fig, ax = plt.subplots(figsize=(10, 4))
    for x, mean, sem, shuffled, first in ((x_l, dg_l, dg_l_sem, dg_sl, True),
                                          (x_r, dg_r, dg_r_sem, dg_sr, False)):
        ax.plot(x, mean, color="steelblue", lw=2, label="Real" if first else None)
        ax.fill_between(x, mean - sem, mean + sem, color="steelblue", alpha=0.12)
        ax.plot(x, shuffled, color="tab:orange", lw=2,
                label="Dinucleotide-shuffled" if first else None)

    ax.axvline(0, color="black", ls="--", lw=1, alpha=0.5, label=f"{mge_name} insertion site")
    ax.set_xlabel(f"Position relative to {mge_name} insertion site")
    ax.set_ylabel("Mean ΔG°₃₇ (kcal/mol)")
    ax.set_title(f"{mge_name}: mean stacking ΔG profile (SantaLucia 1998; "
                 f"n = {len(lefts)} left, {len(rights)} right)")
    ax.annotate('↑ easier to melt', xy=(0.01, 0.95), xycoords='axes fraction',
                fontsize=7, color='grey', va='top')
    ax.legend(fontsize=9)

    fig.tight_layout()
    graphs_dir = Path(output_dir) / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)
    png_path = graphs_dir / f"{mge_name}_MELTING_PROFILE.png"
    fig.savefig(png_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Melting profile plot written to {png_path}")
