"""
DNA 3D curvature analysis using the Bolshoi & Trifonov (1989) wedge model —
the same model underlying DNAcurve. Builds per-sequence 3D backbone coordinates,
computes local curvature, writhe, and bend plane orientation.

Writhe and bend plane are computed on the full insertion site sequence
(left + right flank concatenated) per accession, since these are properties
of the whole target region. Curvature is still plotted positionally.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

logger = logging.getLogger(__name__)

# ----------------------------------------------------------------------------
# Bolshoi & Trifonov (1989) dinucleotide wedge parameters
# Tilt (degrees): rotation about the short axis of the base pair
# Roll (degrees): rotation about the long axis of the base pair
# Twist (degrees): helical twist between successive base pairs
# ----------------------------------------------------------------------------
WEDGE_PARAMS = {
    #          tilt   roll   twist
    "AA": ( 0.00,  0.00, 35.62),
    "AT": ( 0.00, -6.00, 31.50),
    "AC": ( 0.00, -1.50, 34.40),
    "AG": ( 0.00,  3.00, 27.70),
    "TA": ( 0.00,  0.00, 36.00),
    "TT": ( 0.00,  0.00, 35.62),
    "TC": ( 0.00,  6.30, 34.50),
    "TG": ( 0.00,  2.00, 35.00),
    "CA": ( 0.00,  6.30, 34.50),
    "CT": ( 0.00,  3.00, 27.70),
    "CC": ( 0.00,  1.60, 33.67),
    "CG": ( 0.00, -8.00, 29.80),
    "GA": ( 0.00,  3.00, 27.70),
    "GT": ( 0.00, -1.50, 34.40),
    "GC": ( 0.00,  0.00, 36.90),
    "GG": ( 0.00,  1.60, 33.67),
}

RISE_PER_BP = 3.4  # Angstroms


def _rot_x(angle_deg):
    a = np.radians(angle_deg)
    return np.array([
        [1,       0,        0],
        [0,  np.cos(a), -np.sin(a)],
        [0,  np.sin(a),  np.cos(a)],
    ])


def _rot_y(angle_deg):
    a = np.radians(angle_deg)
    return np.array([
        [ np.cos(a), 0, np.sin(a)],
        [       0,   1,       0  ],
        [-np.sin(a), 0, np.cos(a)],
    ])


def _rot_z(angle_deg):
    a = np.radians(angle_deg)
    return np.array([
        [np.cos(a), -np.sin(a), 0],
        [np.sin(a),  np.cos(a), 0],
        [       0,          0,  1],
    ])


def build_backbone(sequence: str) -> np.ndarray:
    """
    Build 3D phosphate backbone coordinates for a DNA sequence
    using the Bolshoi-Trifonov wedge model.

    Returns:
        coords: (N, 3) array of backbone positions in Angstroms, N = len(sequence)
    """
    seq = sequence.upper()
    n = len(seq)
    coords = np.zeros((n, 3))

    position = np.array([0.0, 0.0, 0.0])
    frame = np.eye(3)  # local reference frame (columns = x, y, z axes)
    coords[0] = position

    for i in range(n - 1):
        dinuc = seq[i:i + 2]
        if dinuc not in WEDGE_PARAMS:
            tilt, roll, twist = 0.0, 0.0, 34.3  # neutral fallback for ambiguous bases
        else:
            tilt, roll, twist = WEDGE_PARAMS[dinuc]

        step_rot = _rot_z(twist) @ _rot_y(roll) @ _rot_x(tilt)
        frame = frame @ step_rot
        position = position + frame[:, 2] * RISE_PER_BP
        coords[i + 1] = position

    return coords


def local_curvature(coords: np.ndarray, window: int = 3) -> np.ndarray:
    """
    Compute local curvature at each position as the change in tangent direction
    per unit arc length over a sliding window.

    Returns:
        curvature: array of length N (NaN at boundaries)
    """
    n = len(coords)
    curvature = np.full(n, np.nan)
    half = window // 2

    for i in range(half, n - half):
        v1 = coords[i] - coords[i - half]
        v2 = coords[i + half] - coords[i]
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        if norm1 < 1e-10 or norm2 < 1e-10:
            continue
        v1 /= norm1
        v2 /= norm2
        cross = np.cross(v1, v2)
        sin_angle = np.clip(np.linalg.norm(cross), 0, 1)
        arc = (norm1 + norm2) / 2
        curvature[i] = sin_angle / arc

    return curvature


def writhe(coords: np.ndarray) -> float:
    """
    Compute writhe via the discretised Gauss linking integral.
    Writhe > 0 = right-handed supercoiling, < 0 = left-handed.
    """
    n = len(coords)
    w = 0.0
    for i in range(n - 1):
        for j in range(i + 2, n - 1):
            r1, r2 = coords[i], coords[i + 1]
            r3, r4 = coords[j], coords[j + 1]
            t1 = r2 - r1
            t2 = r4 - r3
            r = (r3 + r4) / 2 - (r1 + r2) / 2
            r_norm = np.linalg.norm(r)
            if r_norm < 1e-10:
                continue
            cross = np.cross(t1, t2)
            w += np.dot(cross, r) / (r_norm ** 3)
    return w / (4 * np.pi)


def bend_plane_angle(coords: np.ndarray) -> float:
    """
    Angle (degrees) between the best-fit plane of the backbone trajectory
    and the XY plane, via PCA. 0° = planar bend in XY, 90° = bend is perpendicular.
    """
    centred = coords - coords.mean(axis=0)
    _, _, vh = np.linalg.svd(centred)
    normal = vh[-1]
    cos_angle = abs(np.dot(normal, [0, 0, 1]))
    return np.degrees(np.arccos(np.clip(cos_angle, 0, 1)))


def _pair_sequences(
    left_seqs: list[str],
    right_seqs: list[str],
    left_names: list[str],
    right_names: list[str],
    flank_len: int,
) -> list[str]:
    """
    Match left and right flanks by accession name and concatenate them
    into full insertion site sequences (left + right, length = 2 * flank_len).

    Accession names are matched after stripping any '_rc' suffix added by
    merge_inverse_symmetric_flanks. Unmatched flanks are silently dropped.

    Returns:
        combined: list of concatenated sequences, one per matched accession pair
    """
    # Strip _rc suffix for matching purposes
    right_map = {
        name.removesuffix("_rc"): seq
        for name, seq in zip(right_names, right_seqs)
        if len(seq) == flank_len
    }

    combined = []
    for name, left_seq in zip(left_names, left_seqs):
        if len(left_seq) != flank_len:
            continue
        clean_name = name.removesuffix("_rc")
        if clean_name in right_map:
            combined.append(left_seq + right_map[clean_name])

    return combined


# ----------------------------------------------------------------------------
# Main analysis function
# ----------------------------------------------------------------------------

def analyse_dna_curvature(
    left_seqs: list[str],
    right_seqs: list[str],
    left_names: list[str],
    right_names: list[str],
    flank_len: int,
    mge_name: str,
    outdir: str,
    curvature_window: int = 5,
):
    """
    Compute and plot DNA curvature metrics at insertion sites.

    Curvature profile: averaged per-position across all sequences, plotted
    with left/right flanks on a common positional axis centred at 0.

    Writhe and bend plane: computed on full insertion site sequences
    (left + right concatenated per accession), since these are whole-molecule
    properties that require the complete target region.

    Outputs (in <outdir>/graphs/):
        curvature_<mge_name>.png   - positional curvature profile
        writhe_<mge_name>.png      - writhe distribution across accessions
        bendplane_<mge_name>.png   - bend plane angle distribution
    """
    graphs_dir = Path(outdir) / "graphs"
    graphs_dir.mkdir(parents=True, exist_ok=True)

    site_len = 2 * flank_len

    # --- Build combined sequences for writhe/bend plane ---
    combined_seqs = _pair_sequences(left_seqs, right_seqs, left_names, right_names, flank_len)
    logger.info(
        f"Paired {len(combined_seqs)}/{len(left_seqs)} left and {len(right_seqs)} right flanks "
        f"into full insertion-site sequences (by accession name)"
    )
    if not combined_seqs:
        logger.warning(
            "No left/right flanks could be paired by accession name -- "
            "writhe and bend-plane results will be empty (curvature profile is unaffected)"
        )
    logger.debug("Building DNA backbone coordinates and computing local curvature...")

    # --- Per-position curvature (still per-flank for positional plot) ---
    def _curvature_matrix(seqs, expected_len):
        matrix = []
        for seq in seqs:
            if len(seq) != expected_len:
                continue
            coords = build_backbone(seq)
            matrix.append(local_curvature(coords, window=curvature_window))
        return np.array(matrix) if matrix else np.full((1, expected_len), np.nan)

    left_curv = _curvature_matrix(left_seqs, flank_len)
    right_curv = _curvature_matrix(right_seqs, flank_len)

    left_mean = np.nanmean(left_curv, axis=0)
    left_sem  = np.nanstd(left_curv, axis=0) / np.sqrt(np.sum(~np.isnan(left_curv), axis=0))
    right_mean = np.nanmean(right_curv, axis=0)
    right_sem  = np.nanstd(right_curv, axis=0) / np.sqrt(np.sum(~np.isnan(right_curv), axis=0))

    left_x  = np.arange(-flank_len, 0)
    right_x = np.arange(1, flank_len + 1)

    # --- Writhe and bend plane on combined sequences ---
    combined_writhes = []
    combined_bends   = []
    for seq in combined_seqs:
        coords = build_backbone(seq)
        combined_writhes.append(writhe(coords))
        combined_bends.append(bend_plane_angle(coords))

    combined_writhes = np.array(combined_writhes)
    combined_bends   = np.array(combined_bends)

    # --- Figure 1: Positional curvature profile ---
    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(left_x, left_mean, color="steelblue", lw=1.8, label="Left flank")
    ax.fill_between(left_x, left_mean - left_sem, left_mean + left_sem,
                    alpha=0.25, color="steelblue")
    ax.plot(right_x, right_mean, color="coral", lw=1.8, label="Right flank")
    ax.fill_between(right_x, right_mean - right_sem, right_mean + right_sem,
                    alpha=0.25, color="coral")
    ax.axvline(0, color="black", linestyle="--", lw=0.9, label=mge_name)
    ax.set_xlabel(f"Position relative to {mge_name}")
    ax.set_ylabel("Mean local curvature (Å⁻¹)")
    ax.set_title(
        f"DNA backbone curvature at {mge_name} insertion sites "
        f"(n={len(left_seqs)} left, {len(right_seqs)} right)"
    )
    ax.legend(framealpha=0.8)
    ax.grid(True, alpha=0.3)
    fig.savefig(graphs_dir / f"{mge_name}_CURVATURE.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Curvature profile plot written to {graphs_dir / f'{mge_name}_CURVATURE.png'}")

    # --- Figure 2: Writhe distribution (combined sequences) ---
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(combined_writhes, bins=30, color="mediumpurple", alpha=0.75)
    ax.axvline(0, color="black", linestyle="--", lw=0.9)
    ax.axvline(np.mean(combined_writhes), color="mediumpurple", linestyle="-",
               lw=1.5, label=f"Mean = {np.mean(combined_writhes):.4f}")
    ax.set_xlabel("Writhe")
    ax.set_ylabel("Count")
    ax.set_title(
        f"Writhe of full insertion site sequences at {mge_name} "
        f"(n={len(combined_seqs)} paired accessions, {site_len} bp each)"
    )
    ax.legend(framealpha=0.8)
    ax.grid(True, alpha=0.3)
    fig.savefig(graphs_dir / f"{mge_name}_WRITHE.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Writhe distribution plot written to {graphs_dir / f'{mge_name}_WRITHE.png'}")

    # --- Figure 3: Bend plane orientation (combined sequences) ---
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(combined_bends, bins=30, color="seagreen", alpha=0.75)
    ax.axvline(np.mean(combined_bends), color="seagreen", linestyle="-",
               lw=1.5, label=f"Mean = {np.mean(combined_bends):.1f}°")
    ax.set_xlabel("Bend plane angle from XY (degrees)")
    ax.set_ylabel("Count")
    ax.set_title(
        f"Bend plane orientation at {mge_name} insertion sites "
        f"(n={len(combined_seqs)} paired accessions)"
    )
    ax.legend(framealpha=0.8)
    ax.grid(True, alpha=0.3)
    fig.savefig(graphs_dir / f"{mge_name}_BENDPLANE.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Bend plane distribution plot written to {graphs_dir / f'{mge_name}_BENDPLANE.png'}")

    result = {
        "n_paired": len(combined_seqs),
        "mean_writhe": float(np.mean(combined_writhes)) if len(combined_writhes) else float("nan"),
        "mean_bend_angle": float(np.mean(combined_bends)) if len(combined_bends) else float("nan"),
    }
    logger.info(
        f"Bendability result: mean writhe = {result['mean_writhe']:.4f}, "
        f"mean bend-plane angle = {result['mean_bend_angle']:.1f}° "
        f"(n={result['n_paired']} paired insertion-site sequences)"
    )
    return result