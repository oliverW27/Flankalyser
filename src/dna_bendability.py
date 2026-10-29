"""
DNA bendability analysis using the Goodsell & Dickerson (1994) BEND model.

Local bend is the angle between base pair normal vectors a fixed number of
steps apart. Normals are built by accumulating the roll/tilt/twist of each
successive step, so the profile reflects how the sequence arranges its
intrinsically bent steps, not just its base composition.

Each real profile is compared against two nulls built from the same sequence:
  - dinucleotide-shuffled -- preserves stacking composition, destroys
    positional arrangement, so a difference here means the *order* matters
  - mononucleotide-shuffled -- preserves base composition only, so a
    difference here means composition alone is doing the work

References:
  Goodsell DS & Dickerson RE (1994) Nucleic Acids Res 22:5497-5503
"""

import logging
import random
import warnings
from contextlib import contextmanager
from functools import lru_cache
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.utils import dinucleotide_shuffle

logger = logging.getLogger(__name__)

MODELS = ["nucleosome", "calladine", "bolshoy", "cacchione",
          "reversed_calladine", "koo_crothers"]

DEFAULT_MODEL = "nucleosome"
DEFAULT_N_B = 1

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


# ----------------------------------------------------------------------------
# Dinucleotide bending models (Table 1 of Goodsell & Dickerson 1994).
# Angles in degrees.
# ----------------------------------------------------------------------------

def _complement_step(step):
    return COMPLEMENT[step[1]] + COMPLEMENT[step[0]]


def _build_dimer_model(roll_tilt_twist):
    """Expand 10 dinucleotide steps into all 16 (complement: roll/twist same, tilt flips)."""
    full = {}
    for step, (roll, tilt, twist) in roll_tilt_twist.items():
        full[step] = (roll, tilt, twist)
        rc = _complement_step(step)
        if rc not in full:
            full[rc] = (roll, -tilt, twist)
    return full


_CALLADINE = _build_dimer_model({
    "AA": (0.0, 0.0, 35.0), "GA": (3.3, 0.0, 34.0),
    "AG": (3.3, 0.0, 34.0), "GG": (3.3, 0.0, 34.0),
    "AT": (3.3, 0.0, 34.0), "AC": (3.3, 0.0, 34.0),
    "GC": (3.3, 0.0, 34.0), "TA": (6.6, 0.0, 34.0),
    "CA": (3.3, 0.0, 34.0), "CG": (3.3, 0.0, 34.0),
})

_BOLSHOY = _build_dimer_model({
    "AA": (-6.5,  3.2, 35.62), "GA": (-2.7, -4.6, 36.90),
    "AG": ( 8.4, -0.3, 27.70), "GG": ( 1.2, -1.8, 33.67),
    "AT": ( 2.6,  0.0, 31.50), "AC": (-0.9, -0.7, 34.40),
    "GC": (-5.0,  0.0, 40.00), "TA": ( 0.9,  0.0, 36.00),
    "CA": ( 1.6,  3.1, 34.50), "CG": ( 6.7,  0.0, 29.80),
})

_CACCHIONE = _build_dimer_model({
    "AA": (-5.4, -0.5, 35.9), "GA": ( 2.0, -1.7, 35.8),
    "AG": ( 1.0, -1.6, 35.6), "GG": ( 1.3, -0.6, 33.0),
    "AT": (-7.3,  0.0, 35.0), "AC": (-2.4, -2.7, 34.6),
    "GC": (-3.7,  0.0, 33.3), "TA": ( 8.0,  0.0, 34.6),
    "CA": ( 6.7,  0.4, 34.5), "CG": ( 4.6,  0.0, 33.7),
})

_REVERSED_CALLADINE = _build_dimer_model({
    "AA": ( 3.3, 0.0, 35.0), "GA": (0.0, 0.0, 34.0),
    "AG": ( 0.0, 0.0, 34.0), "GG": (0.0, 0.0, 34.0),
    "AT": ( 0.0, 0.0, 34.0), "AC": (0.0, 0.0, 34.0),
    "GC": ( 0.0, 0.0, 34.0), "TA": (-3.3, 0.0, 34.0),
    "CA": ( 0.0, 0.0, 34.0), "CG": (0.0, 0.0, 34.0),
})

_DIMER_TABLES = {
    "calladine": _CALLADINE,
    "bolshoy": _BOLSHOY,
    "cacchione": _CACCHIONE,
    "reversed_calladine": _REVERSED_CALLADINE,
}


# ----------------------------------------------------------------------------
# Nucleosome positioning model (trimer)
# ----------------------------------------------------------------------------

# Percentage of each triplet facing outwards from the nucleosome core, i.e.
# how strongly it prefers to roll away from the histone surface.
_PERCENT_OUT = {
    "GGC": 45, "TCG": 31, "AGC": 25, "CGC": 25, "CAT": 18, "CAC": 17,
    "GGG": 13, "TGC": 13, "AGT": 11, "GAG":  8, "GGT":  8, "TGG":  8,
    "CGT":  8, "AGG":  8, "GAC":  8, "TGA":  8, "GAT":  7, "TGT":  6,
    "AAG":  6, "CGG":  2, "CAG": -2, "GGA": -5, "TAC": -6, "AAC": -6,
    "AGA": -9, "CAA": -9, "GAA": -12, "TAT": -13, "TAG": -18, "TAA": -20,
    "AAT": -30, "AAA": -36,
}


def _trimer_roll(triplet, max_roll=10.0):
    """Convert percent-out to roll using the linear scaling from the paper."""
    rc = "".join(COMPLEMENT[b] for b in triplet[::-1])
    pct = _PERCENT_OUT.get(triplet, _PERCENT_OUT.get(rc))
    if pct is None:
        return 0.0
    return max_roll * (pct + 36) / (45 + 36)


# ----------------------------------------------------------------------------
# Koo & Crothers junction model
# ----------------------------------------------------------------------------

def _find_runs(seq, base, min_len=3):
    """Return inclusive (start, end) indices of runs of `base` of length >= min_len."""
    out = []
    n = len(seq)
    j = 0
    while j < n:
        if seq[j] == base:
            k = j
            while k < n and seq[k] == base:
                k += 1
            if k - j >= min_len:
                out.append((j, k - 1))
            j = k
        else:
            j += 1
    return out


def _koo_crothers_step_params(seq, i):
    """Roll/tilt/twist for step seq[i] -> seq[i+1] under the Koo & Crothers model.

    Bending is localised to the junctions where an A-tract (or its T-tract
    mirror) meets the surrounding sequence, rather than spread over every step.
    """
    twist = 34.3

    for start, end in _find_runs(seq, "A"):
        if i == start - 1 and start > 0:
            prev = seq[start - 1]
            if prev in ("G", "C"):
                return (5.8, -9.6, twist)
            if prev == "T":
                return (11.6, 0.0, twist)
        if i == end and end < len(seq) - 1:
            nxt = seq[end + 1]
            if nxt in ("G", "C"):
                return (0.0, 12.0, twist)
            if nxt == "T":
                return (0.0, 0.0, twist)

    for start, end in _find_runs(seq, "T"):
        if i == end and end < len(seq) - 1:
            nxt = seq[end + 1]
            if nxt in ("C", "G"):
                return (5.8, 9.6, twist)
            if nxt == "A":
                return (11.6, 0.0, twist)
        if i == start - 1 and start > 0:
            prev = seq[start - 1]
            if prev in ("C", "G"):
                return (0.0, -12.0, twist)
            if prev == "A":
                return (0.0, 0.0, twist)

    return (0.0, 0.0, twist)


# ----------------------------------------------------------------------------
# Rotation helpers
# ----------------------------------------------------------------------------

def _rot_x(deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rot_y(deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def _rot_z(deg):
    a = np.radians(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


@lru_cache(maxsize=None)
def _step_rotation(roll, tilt, twist):
    """Rotation matrix for one base pair step, cached -- the models only ever
    emit a handful of distinct (roll, tilt, twist) triples, and this runs
    inside a per-base loop over thousands of sequences."""
    return _rot_z(twist) @ _rot_y(tilt) @ _rot_x(roll)


def _step_params(seq, i, model):
    """Return (roll, tilt, twist) in degrees for the step seq[i] -> seq[i+1]."""
    if model == "nucleosome":
        triplet = seq[i:i + 3] if i + 2 < len(seq) else seq[i - 1:i + 2]
        return (_trimer_roll(triplet), 0.0, 34.3)
    if model == "koo_crothers":
        return _koo_crothers_step_params(seq, i)
    if model in _DIMER_TABLES:
        return _DIMER_TABLES[model].get(seq[i:i + 2], (0.0, 0.0, 34.0))
    raise ValueError(f"Unknown bending model: {model!r}. Available: {MODELS}")


# ----------------------------------------------------------------------------
# Bend profile
# ----------------------------------------------------------------------------

def compute_normals(seq, model=DEFAULT_MODEL):
    """Build the local frame at each base pair by accumulating twist/roll/tilt.

    Returns:
        normals: (N, 3) array of base pair normal vectors, N = len(seq)
    """
    seq = seq.upper()
    n = len(seq)
    frames = np.zeros((n, 3, 3))
    frames[0] = np.eye(3)

    for i in range(n - 1):
        frames[i + 1] = frames[i] @ _step_rotation(*_step_params(seq, i, model))

    return frames[:, :, 2]


def bend_profile(seq, model=DEFAULT_MODEL, n_b=DEFAULT_N_B):
    """Local bend angle (degrees) at each position of a sequence.

    Args:
        seq: DNA sequence (ACGT, case-insensitive)
        model: one of MODELS
        n_b: bracket half-width -- bend at position i is the angle between the
             normals at i - n_b and i + n_b (the paper recommends 1 or 2)

    Returns:
        bend: array of length len(seq), NaN in the first and last n_b positions
              where the bracket runs off the end
    """
    seq = seq.upper()
    n = len(seq)
    normals = compute_normals(seq, model)

    bend = np.full(n, np.nan)
    for i in range(n_b, n - n_b):
        cos_angle = np.clip(np.dot(normals[i - n_b], normals[i + n_b]), -1.0, 1.0)
        bend[i] = np.degrees(np.arccos(cos_angle))

    return bend


@contextmanager
def _quiet_empty_slices():
    """Suppress the all-NaN-slice warnings from column-wise nanmean/nanstd.

    The first and last n_b positions of every profile are NaN by construction,
    so those columns are legitimately empty -- the NaN that comes back is the
    right answer, not something to warn about.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Mean of empty slice", RuntimeWarning)
        warnings.filterwarnings("ignore", "Degrees of freedom <= 0", RuntimeWarning)
        yield


def mono_shuffle(seq):
    """Shuffle a sequence preserving base composition only."""
    bases = list(seq)
    random.shuffle(bases)
    return "".join(bases)


# ----------------------------------------------------------------------------
# Input pairing
# ----------------------------------------------------------------------------

def pair_flanks(left_data, right_data, acc_names, flank_len):
    """Pair left and right flanks from the same accession into insertion-site
    sequences (left + right, length 2 * flank_len).

    Unlike the per-side deduplication done elsewhere in the pipeline, pairing
    has to happen before deduplication -- bend is a property of the joined
    sequence, so the left/right correspondence on each input line must be kept.
    Pipe-separated segments are zipped positionally, matching how TETyper emits
    alternative flank calls for the same site.

    Returns:
        pairs: list of (accession, combined_sequence), deduplicated on sequence
    """
    seen, pairs, skipped = set(), [], 0

    for name, left_field, right_field in zip(acc_names, left_data, right_data):
        for left, right in zip(left_field.upper().split("|"),
                               right_field.upper().split("|")):
            left, right = left.strip(), right.strip()
            if len(left) != flank_len or len(right) != flank_len:
                skipped += 1
                continue
            combined = left + right
            if any(base not in "ACGT" for base in combined):
                skipped += 1
                continue
            if combined not in seen:
                seen.add(combined)
                pairs.append((name, combined))

    logger.info(
        f"Paired {len(pairs)} unique insertion-site sequences of {2 * flank_len} bp "
        f"from {len(acc_names)} accessions ({skipped} segment(s) skipped: wrong "
        f"length or non-ACGT)"
    )
    return pairs


# ----------------------------------------------------------------------------
# Main analysis function
# ----------------------------------------------------------------------------

def run_bendability_analysis(
    left_data: list[str],
    right_data: list[str],
    acc_names: list[str],
    flank_len: int,
    mge_name: str,
    output_dir: str,
    model: str = DEFAULT_MODEL,
    n_b: int = DEFAULT_N_B,
    n_shuffles: int = 5,
    max_seqs: int = 5000,
    seed: int | None = None,
):
    """Compute and plot the mean local bend profile across insertion sites.

    Left and right flanks from each accession are joined into a single
    insertion-site sequence, so the midpoint of the profile is the insertion
    point itself. The real profile is plotted against dinucleotide- and
    mononucleotide-shuffled nulls generated from the same sequences, each as
    mean +/- SD across sequences.

    Args:
        left_data, right_data, acc_names: raw per-line fields from the summary
            file (pipe-separated segments are handled here)
        flank_len: length of each flank, so each joined sequence is 2 * flank_len
        mge_name: name of the MGE, used in the title and output filenames
        output_dir: writes <output_dir>/graphs/<mge_name>_BEND_PROFILE.png and
            <output_dir>/files/<mge_name>_bend_profile.tsv
        model: bending model, one of MODELS
        n_b: bracket half-width for the local bend angle
        n_shuffles: shuffles per sequence per null, averaged before pooling
        max_seqs: downsample to at most this many sequences (0 for no limit)
        seed: random seed for downsampling and shuffling, for reproducible runs

    Returns:
        dict with n_seqs and the mean real/dinuc/mono bend, or None if no
        usable sequences were found.
    """
    if model not in MODELS:
        raise ValueError(f"Unknown bending model: {model!r}. Available: {MODELS}")

    if seed is not None:
        random.seed(seed)

    pairs = pair_flanks(left_data, right_data, acc_names, flank_len)
    if not pairs:
        logger.warning(
            "No left/right flank pairs survived filtering -- skipping bendability analysis"
        )
        return None

    if max_seqs and len(pairs) > max_seqs:
        pairs = random.sample(pairs, max_seqs)
        logger.info(f"Downsampled to {max_seqs} sequences for the bend profile")

    logger.info(
        f"Computing local bend profiles (model={model}, n_b={n_b}) with "
        f"{n_shuffles} dinucleotide + {n_shuffles} mononucleotide shuffle(s) per sequence..."
    )

    real_profiles, dinuc_profiles, mono_profiles = [], [], []
    n_shuffle_failed = 0

    for _name, seq in pairs:
        dinuc_seqs = []
        for _ in range(n_shuffles):
            try:
                dinuc_seqs.append(dinucleotide_shuffle(seq))
            except RuntimeError:
                continue
        if not dinuc_seqs:
            n_shuffle_failed += 1
            continue

        real_profiles.append(bend_profile(seq, model, n_b))
        with _quiet_empty_slices():
            dinuc_profiles.append(np.nanmean(
                [bend_profile(s, model, n_b) for s in dinuc_seqs], axis=0))
            mono_profiles.append(np.nanmean(
                [bend_profile(mono_shuffle(seq), model, n_b) for _ in range(n_shuffles)],
                axis=0))

    if n_shuffle_failed:
        logger.warning(
            f"{n_shuffle_failed} sequence(s) dropped -- dinucleotide shuffle failed"
        )
    if not real_profiles:
        logger.warning("No bend profiles could be computed -- skipping bendability analysis")
        return None

    real_arr = np.array(real_profiles)
    dinuc_arr = np.array(dinuc_profiles)
    mono_arr = np.array(mono_profiles)

    graphs_dir = Path(output_dir) / "graphs"
    files_dir = Path(output_dir) / "files"
    graphs_dir.mkdir(parents=True, exist_ok=True)
    files_dir.mkdir(parents=True, exist_ok=True)

    means = _plot_mean_bend(real_arr, dinuc_arr, mono_arr, model, mge_name, graphs_dir)
    _write_profile_tsv(means, mge_name, files_dir)

    with _quiet_empty_slices():
        result = {
            "n_seqs": len(real_arr),
            "mean_real_bend": float(np.nanmean(real_arr)),
            "mean_dinuc_bend": float(np.nanmean(dinuc_arr)),
            "mean_mono_bend": float(np.nanmean(mono_arr)),
        }
    logger.info(
        f"Bendability result: mean local bend = {result['mean_real_bend']:.3f} deg "
        f"(dinuc-shuffled {result['mean_dinuc_bend']:.3f}, "
        f"mono-shuffled {result['mean_mono_bend']:.3f}) "
        f"across n={result['n_seqs']} insertion-site sequences"
    )
    return result


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------

def _plot_mean_bend(real_arr, dinuc_arr, mono_arr, model, mge_name, graphs_dir):
    """Plot mean +/- SD bend for the real sequences and both shuffled nulls.

    Returns the per-position (mean, sd) of each series, for the TSV.
    """
    n = len(real_arr)
    x = np.arange(real_arr.shape[1])

    series = [
        ("real", real_arr, "steelblue", f"Real (n={n})"),
        ("dinuc", dinuc_arr, "mediumseagreen", f"Dinuc-shuffled (n={n})"),
        ("mono", mono_arr, "darkorange", f"Mono-shuffled (n={n})"),
    ]

    fig, ax = plt.subplots(figsize=(10, 4))
    means = {}
    for key, arr, colour, label in series:
        with _quiet_empty_slices():
            mean = np.nanmean(arr, axis=0)
            sd = np.nanstd(arr, axis=0)
        means[key] = (mean, sd)
        ax.fill_between(x, mean - sd, mean + sd, alpha=0.2, color=colour)
        ax.plot(x, mean, color=colour, lw=2, label=label)

    ax.axvline(real_arr.shape[1] / 2.0, color="black", ls="--", lw=1, alpha=0.5,
               label=f"{mge_name} insertion site")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel(f"Local bend (degrees, {model} model)")
    ax.set_title(f"{mge_name}: mean bend profile")
    ax.legend(fontsize=9)
    fig.tight_layout()

    png_path = graphs_dir / f"{mge_name}_BEND_PROFILE.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Mean bend profile plot written to {png_path}")
    return means


def _write_profile_tsv(means, mge_name, files_dir):
    """Write the plotted per-position means and SDs so the profile can be
    re-plotted or compared across elements without recomputing."""
    tsv_path = files_dir / f"{mge_name}_bend_profile.tsv"
    with open(tsv_path, "w") as f:
        f.write("position\treal_mean\treal_sd\tdinuc_mean\tdinuc_sd\t"
                "mono_mean\tmono_sd\n")
        for i in range(len(means["real"][0])):
            f.write(f"{i}\t"
                    f"{means['real'][0][i]:.4f}\t{means['real'][1][i]:.4f}\t"
                    f"{means['dinuc'][0][i]:.4f}\t{means['dinuc'][1][i]:.4f}\t"
                    f"{means['mono'][0][i]:.4f}\t{means['mono'][1][i]:.4f}\n")
    logger.info(f"Bend profile table written to {tsv_path}")
