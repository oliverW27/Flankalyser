"""
Nucleotide Distribution Analysis Module
Analyzes position-specific nucleotide frequencies and statistical significance.
"""

import matplotlib.pyplot as plt
from scipy.stats import chisquare
import numpy as np
import logging

logger = logging.getLogger(__name__)


def count_nucl(l_data, r_data, FLANK_LEN):
    """Count nucleotides at each position for left and right flanks."""
    logger.debug("Counting nucleotides for flank sequences...")
    l_counts = {"A": [0]*FLANK_LEN, "C": [0]*FLANK_LEN, "G": [0]*FLANK_LEN, "T": [0]*FLANK_LEN}
    r_counts = {"A": [0]*FLANK_LEN, "C": [0]*FLANK_LEN, "G": [0]*FLANK_LEN, "T": [0]*FLANK_LEN}
    
    for seq in l_data:
        if len(seq) != FLANK_LEN:
            logger.warning(f"Left flank sequence length mismatch: expected {FLANK_LEN}, got {len(seq)}")
        for i, nucl in enumerate(seq):
            if nucl in l_counts:
                l_counts[nucl][i] += 1
            elif nucl != "N":
                logger.warning(f"Unexpected nucleotide '{nucl}' in left flank at position {i}")
                
    for seq in r_data:
        if len(seq) != FLANK_LEN:
            logger.warning(f"Right flank sequence length mismatch: expected {FLANK_LEN}, got {len(seq)}")
        for i, nucl in enumerate(seq):
            if nucl in r_counts:
                r_counts[nucl][i] += 1
            elif nucl != "N":
                logger.warning(f"Unexpected nucleotide '{nucl}' in right flank at position {i}")
                
    logger.debug("Nucleotide counting complete.")
    return l_counts, r_counts


def nucl_dist_stat_sig(data):
    """Calculate chi-square test p-values for nucleotide distribution at each position.

    Args:
        data: Dict with counts for A, C, G, T at each position

    Returns:
        List of p-values for each position
    """
    logger.debug("Calculating chi-square test p-values...")

    expected_props = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    observed = np.array([data['A'], data['C'], data['G'], data['T']])
    total_counts = np.sum(observed, axis=0)
    
    # Calculate expected counts using GC-adjusted proportions
    expected = np.zeros_like(observed, dtype=float)
    for i, nucl in enumerate(['A', 'C', 'G', 'T']):
        expected[i] = total_counts * expected_props[nucl]
    
    chi2_stats = []
    p_values = []
    
    for i in range(observed.shape[1]):
        obs = observed[:, i]
        exp = expected[:, i]
        try:
            chi2_stat, p_value = chisquare(f_obs=obs, f_exp=exp)
            chi2_stats.append(chi2_stat)
            p_values.append(p_value)
        except ValueError as e:
            logger.error(f"Chi-square test failed at position {i}: {e}")
            chi2_stats.append(np.nan)
            p_values.append(np.nan)
    
    logger.debug("Chi-square calculations complete.")
    return p_values


def calc_nucl_props(l_data, r_data, FLANK_LEN):
    """Calculate nucleotide proportions combining left and right flank data."""
    logger.debug("Calculating combined nucleotide proportions...")
    combined_dict = {}
    for nucleotide in l_data:
        combined_dict[nucleotide] = l_data[nucleotide] + [np.nan] + r_data[nucleotide]

    combined_proportions = {nucleotide: [] for nucleotide in combined_dict}
    num_loci = len(combined_dict['A'])

    for i in range(num_loci):
        total_count = sum(combined_dict[n][i] for n in combined_dict)
        for n in combined_dict:
            if total_count != 0:
                proportion = combined_dict[n][i] / total_count
            else:
                proportion = np.nan
                logger.debug(f"No counts at locus {i}, assigning NaN proportion")
            combined_proportions[n].append(proportion)
            
    logger.debug("Nucleotide proportion calculation complete.")
    return combined_proportions



def get_stats_ast(l_stats, r_stats, n):
    """Generate significance asterisks based on p-values with Bonferroni correction."""
    logger.debug("Generating significance asterisks with Bonferroni correction...")
    ast = []
    for i in l_stats + ["NaN"] + r_stats:
        if i == "NaN":
            ast.append(" ")
        elif i < 0.001/n:
            ast.append("***")
        elif i < 0.01/n:
            ast.append("**")
        elif i < 0.05/n:
            ast.append("*")
        else:
            ast.append(" ")
    return ast


def graph(
    FLANK_LEN,
    IS_NAME,
    combined_proportions,
    l_stats,
    r_stats,
    outdir,
    total_left=None,
    unique_left=None,
    total_right=None,
    unique_right=None,
):
    """Generate nucleotide frequency plot with significance indicators and flank counts."""
    logger.debug(f"Generating nucleotide frequency graph for {IS_NAME}...")
    stats_ast = get_stats_ast(l_stats, r_stats, 2*FLANK_LEN)
    x = list(range(-FLANK_LEN, FLANK_LEN + 1))
    colours = {"A": "green", "C": "blue", "G": "orange", "T": "red"}

    for nt in combined_proportions:
        combined_proportions[nt] = [
            np.nan if str(v).lower() == 'nan' else v
            for v in combined_proportions[nt]
        ]

    plt.figure(figsize=(24, 14))

    for nucleotide in combined_proportions:
        plt.plot(
            x,
            combined_proportions[nucleotide],
            label=nucleotide,
            color=colours[nucleotide],
            linewidth=5,
        )

    # Add horizontal line for expected proportion (25% each)
    plt.axhline(y=0.25, color='gray', linestyle=':', linewidth=3,
                label='Expected (25%)')

    max_y_value = np.nanmax([
        np.nanmax(combined_proportions["A"]),
        np.nanmax(combined_proportions["C"]),
        np.nanmax(combined_proportions["G"]),
        np.nanmax(combined_proportions["T"]),
    ])
    plt.ylim(0, max_y_value + 0.12)

    for i, xpos in enumerate(x):
        max_nuc_value = max(
            combined_proportions["A"][i],
            combined_proportions["C"][i],
            combined_proportions["G"][i],
            combined_proportions["T"][i],
        )
        y_ast = max_nuc_value + 0.01
        ast = stats_ast[i]
        if ast != " ":
            for _ in range(len(ast)):
                plt.text(xpos, y_ast, "*", ha="center", va="bottom", fontsize=18)
                y_ast += 0.005

    plt.axvline(x=0, color="gray", linestyle="--", label=IS_NAME, linewidth=1.5)

    plt.xlabel(f"Position relative to {IS_NAME}", fontsize=25)
    plt.ylabel("Proportional frequency", fontsize=25)
    plt.title(f"Frequency of nucleotides in {IS_NAME} flanking sequences ({total_left}, {unique_left} | {total_right}, {unique_right})", fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(fontsize=25)
    plt.grid(True)
    plt.tight_layout()

    outpath = f"{outdir}/graphs/{IS_NAME}_DNA_FREQ.png"
    try:
        plt.savefig(outpath)
        logger.info(f"Saved graph to {outpath}")
    except Exception as e:
        logger.exception(f"Failed to save graph: {e}")