"""
DNA bendability analysis using the nucleosome positioning (trimer) model
from Goodsell & Dickerson (1994).

Local bend at a position is the angle between the base pair normals one step
either side of it. Normals are built by accumulating the roll and twist of
each successive step, so the profile reflects how the sequence arranges its
intrinsically bent steps, not just its base composition.

Left and right flanks are profiled separately, each against a
dinucleotide-shuffled control (same stacking composition, order destroyed),
so a difference between the two lines means the order of the sequence matters.
The flanks are not joined, because the insertion (e.g. a target site
duplication) means the joined sequence is not the original target site.

References:
  Goodsell DS & Dickerson RE (1994) Nucleic Acids Res 22:5497-5503
"""

import logging
import random
from functools import lru_cache
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.utils import dinucleotide_shuffle, reverse_complement

logger = logging.getLogger(__name__)

TWIST = 34.3      # degrees per step, constant in this model
MAX_ROLL = 10.0   # roll (degrees) given to the most outward-facing trinucleotide
MAX_SEQS = 5000   # downsample each flank to this many sequences to keep runtime reasonable

# Percentage of each triplet facing outwards from the nucleosome core, i.e.
# how strongly it prefers to roll away from the histone surface. Only one of
# each reverse-complement pair is listed.
_PERCENT_OUT = {
    "GGC": 45, "TCG": 31, "AGC": 25, "CGC": 25, "CAT": 18, "CAC": 17,
    "GGG": 13, "TGC": 13, "AGT": 11, "GAG":  8, "GGT":  8, "TGG":  8,
    "CGT":  8, "AGG":  8, "GAC":  8, "TGA":  8, "GAT":  7, "TGT":  6,
    "AAG":  6, "CGG":  2, "CAG": -2, "GGA": -5, "TAC": -6, "AAC": -6,
    "AGA": -9, "CAA": -9, "GAA": -12, "TAT": -13, "TAG": -18, "TAA": -20,
    "AAT": -30, "AAA": -36,
}


def _trimer_roll(triplet):
    """Convert percent-out to roll, scaled linearly from 0 (AAA) to MAX_ROLL (GGC)."""
    pct = _PERCENT_OUT.get(triplet, _PERCENT_OUT.get(reverse_complement(triplet)))
    return MAX_ROLL * (pct + 36) / (45 + 36)


@lru_cache(maxsize=None)
def _step_rotation(roll):
    """Rotation for one base pair step: roll about x, then twist about z.

    Cached because there are only 21 distinct roll values, and this runs
    inside a per-base loop over thousands of sequences.
    """
    r, t = np.radians(roll), np.radians(TWIST)
    roll_m = np.array([[1, 0, 0], [0, np.cos(r), -np.sin(r)], [0, np.sin(r), np.cos(r)]])
    twist_m = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    return twist_m @ roll_m


def bend_profile(seq):
    """Local bend angle (degrees) at positions 1 .. len(seq) - 2 of an ACGT sequence."""
    frame = np.eye(3)
    normals = [frame[:, 2]]
    for i in range(len(seq) - 1):
        # The last step has no trinucleotide ahead of it, so reuse the one ending at i + 1
        triplet = seq[i:i + 3] if i + 2 < len(seq) else seq[i - 1:i + 2]
        frame = frame @ _step_rotation(_trimer_roll(triplet))
        normals.append(frame[:, 2])
    normals = np.array(normals)
    cos_angle = (normals[:-2] * normals[2:]).sum(axis=1)
    return np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))


def _mean_profiles(seqs, flank_len):
    """Mean bend profile of seqs and of one dinucleotide shuffle of each.
    Returns (real_mean, shuffled_mean, n), or None if nothing could be shuffled."""
    real_sum = np.zeros(flank_len - 2)
    shuf_sum = np.zeros(flank_len - 2)
    n = 0
    for seq in seqs:
        try:
            shuffled = dinucleotide_shuffle(seq)
        except RuntimeError:
            continue
        real_sum += bend_profile(seq)
        shuf_sum += bend_profile(shuffled)
        n += 1
    return (real_sum / n, shuf_sum / n, n) if n else None


def run_bendability_analysis(left_seqs, right_seqs, flank_len, mge_name, output_dir):
    """Plot the mean local bend profile of each flank against a dinucleotide-shuffled
    control, to <output_dir>/graphs/<mge_name>_BEND_PROFILE.png.

    left_seqs / right_seqs are the unique flank sequences from flankalyser.py's
    extract_unique_sequences. Flanks containing anything other than A/C/G/T are skipped.
    """
    fig, ax = plt.subplots(figsize=(10, 4))
    # Bend at a base needs the base either side, so each flank loses its outermost base
    # at both ends: left covers -(flank_len - 1)..-2, right covers +2..+(flank_len - 1)
    sides = (("left", left_seqs, np.arange(1, flank_len - 1) - flank_len),
             ("right", right_seqs, np.arange(1, flank_len - 1) + 1))
    n_used = {}
    for side, seqs, x in sides:
        seqs = [s.upper() for s in seqs if set(s.upper()) <= set("ACGT")]
        if len(seqs) > MAX_SEQS:
            seqs = random.sample(seqs, MAX_SEQS)
            logger.info(f"Downsampled {side} flank to {MAX_SEQS} sequences for the bend profile")
        result = _mean_profiles(seqs, flank_len) if seqs else None
        if result is None:
            logger.warning(f"No usable {side} flank sequences -- leaving it off the bend profile")
            continue
        real, shuffled, n_used[side] = result
        first = len(n_used) == 1  # only label the first side plotted, so the legend has one entry each
        ax.plot(x, real, color="steelblue", lw=2, label="Real" if first else None)
        ax.plot(x, shuffled, color="tab:orange", lw=2,
                label="Dinucleotide-shuffled" if first else None)

    if not n_used:
        logger.warning("No usable flank sequences on either side -- skipping bendability analysis")
        plt.close(fig)
        return

    ax.axvline(0, color="black", ls="--", lw=1, alpha=0.5, label=f"{mge_name} insertion site")
    ax.set_xlabel(f"Position relative to {mge_name} insertion site")
    ax.set_ylabel("Mean local bend (degrees)")
    counts = ", ".join(f"{n} {side}" for side, n in n_used.items())
    ax.set_title(f"{mge_name}: mean bend profile (nucleosome model; n = {counts})")
    ax.legend(fontsize=9)
    fig.tight_layout()

    graphs_dir = Path(output_dir) / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)
    png_path = graphs_dir / f"{mge_name}_BEND_PROFILE.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Mean bend profile plot written to {png_path}")
