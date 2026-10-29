"""
Nucleotide Distribution Analysis Module
Analyses position-specific nucleotide frequencies and statistical significance.
"""

import matplotlib.pyplot as plt
from scipy.stats import chisquare
import numpy as np
import logging

logger = logging.getLogger(__name__)

BASES = "ACGT"


def _count_side(seqs, flank_len, side):
    """Per-position A/C/G/T counts for one flank. Ns are skipped silently."""
    counts = {b: [0] * flank_len for b in BASES}
    for seq in seqs:
        for i, nucl in enumerate(seq):
            if nucl in counts:
                counts[nucl][i] += 1
            elif nucl != "N":
                logger.warning(f"Unexpected nucleotide '{nucl}' in {side} flank at position {i}")
    return counts


def count_nucl(l_data, r_data, flank_len):
    """Count nucleotides at each position for left and right flanks."""
    return _count_side(l_data, flank_len, "left"), _count_side(r_data, flank_len, "right")


def expected_proportions(gc_content=None):
    """Expected A/C/G/T proportions: flat 25% each, or split by GC content."""
    if gc_content is None:
        gc_content = 0.5  # i.e. 25% for each base
    at, gc = (1 - gc_content) / 2, gc_content / 2
    return {'A': at, 'C': gc, 'G': gc, 'T': at}


def nucl_dist_stat_sig(data, gc_content=None):
    """Calculate chi-square test p-values for nucleotide distribution at each position.

    Args:
        data: Dict with counts for A, C, G, T at each position
        gc_content: Expected GC content as decimal (e.g., 0.52 for 52%), used to set
                    the chi-square test's expected proportions to a given background
                    distribution instead of a flat 25% each. If None, assumes 0.25.

    Returns:
        List of p-values for each position
    """
    if gc_content is not None:
        logger.debug(f"Using GC-adjusted expectations: A/T={(1 - gc_content) / 2:.3f}, "
                     f"C/G={gc_content / 2:.3f}")
    expected_props = np.array([expected_proportions(gc_content)[b] for b in BASES])

    observed = np.array([data[b] for b in BASES])
    expected = expected_props[:, None] * observed.sum(axis=0)
    return list(chisquare(f_obs=observed, f_exp=expected, axis=0).pvalue)


def calc_nucl_props(l_data, r_data):
    """Per-position nucleotide proportions across left flank, a NaN gap for the
    insertion site, then right flank. Positions with no counts are NaN."""
    counts = np.array([l_data[b] + [np.nan] + r_data[b] for b in BASES], dtype=float)
    with np.errstate(invalid="ignore"):
        props = counts / counts.sum(axis=0)
    return dict(zip(BASES, props))


def get_stats_ast(l_stats, r_stats, n):
    """Generate significance asterisks based on p-values with Bonferroni correction.
    The insertion-site gap (and any NaN p-value) gets no asterisks."""
    ast = []
    for p in l_stats + [np.nan] + r_stats:
        if p < 0.001/n:
            ast.append("***")
        elif p < 0.01/n:
            ast.append("**")
        elif p < 0.05/n:
            ast.append("*")
        else:
            ast.append(" ")
    return ast


def graph_nucl_dist(
    flank_len,
    mge_name,
    combined_proportions,
    l_stats,
    r_stats,
    outdir,
    total_left,
    unique_left,
    total_right,
    unique_right,
    gc_content=None,
):
    """Generate nucleotide frequency plot with significance indicators and flank counts."""
    stats_ast = get_stats_ast(l_stats, r_stats, 2*flank_len)
    x = list(range(-flank_len, flank_len + 1))
    colours = {"A": "green", "C": "blue", "G": "orange", "T": "red"}

    plt.figure(figsize=(24, 14))

    for nucleotide in combined_proportions:
        plt.plot(
            x,
            combined_proportions[nucleotide],
            label=nucleotide,
            color=colours[nucleotide],
            linewidth=5,
        )

    # Add horizontal lines for the expected background proportions
    if gc_content is not None:
        gc_half = gc_content / 2
        at_half = (1 - gc_content) / 2
        plt.axhline(y=gc_half, color='orange', linestyle='--', linewidth=3,
                    label=f'Expected G/C ({gc_half:.1%})')
        plt.axhline(y=at_half, color='green', linestyle=':', linewidth=3,
                    label=f'Expected A/T ({at_half:.1%})')
    else:
        plt.axhline(y=0.25, color='gray', linestyle=':', linewidth=3,
                    label='Expected (25%)')

    props = np.array([combined_proportions[b] for b in BASES])
    plt.ylim(0, np.nanmax(props) + 0.12)

    # Stack asterisks just above the highest base at each significant position
    for i, xpos in enumerate(x):
        y_ast = props[:, i].max() + 0.01
        if stats_ast[i] != " ":
            for _ in range(len(stats_ast[i])):
                plt.text(xpos, y_ast, "*", ha="center", va="bottom", fontsize=18)
                y_ast += 0.005

    plt.axvline(x=0, color="gray", linestyle="--", label=mge_name, linewidth=1.5)

    plt.xlabel(f"Position relative to {mge_name}", fontsize=25)
    plt.ylabel("Proportional frequency", fontsize=25)
    plt.title(f"Frequency of nucleotides in {mge_name} flanking sequences ({total_left}, {unique_left} | {total_right}, {unique_right})", fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=25)
    plt.grid(True)
    plt.tight_layout()

    outpath = f"{outdir}/graphs/{mge_name}_DNA_FREQ.png"
    try:
        plt.savefig(outpath)
        logger.info(f"Saved graph to {outpath}")
    except Exception as e:
        logger.exception(f"Failed to save graph: {e}")
    plt.close()
