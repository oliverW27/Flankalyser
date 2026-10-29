"""
Sequence logos for the left and right flanks, using logomaker.
This is somewhat redundant with the nucleotide distribution graph, but it makes
the consensus sequence easier to read.
"""
import logging

import pandas as pd
import matplotlib.pyplot as plt
import logomaker

logger = logging.getLogger(__name__)


def generate_sequence_logo(left_counts, right_counts, output_dir, mge_name):
    """Plot a sequence logo for each flank.

    Args:
        left_counts, right_counts: per-position base counts from nucl_dist.count_nucl,
            i.e. {'A': [...], 'C': [...], 'G': [...], 'T': [...]}
        output_dir: writes <output_dir>/graphs/<mge_name>_{LEFT,RIGHT}_LOGO.png
        mge_name: name of the MGE, used as the output filename prefix
    """
    for side, counts in (("LEFT", left_counts), ("RIGHT", right_counts)):
        counts = pd.DataFrame(counts)  # rows = positions, columns = A, C, G, T

        # A position that is N in every flank has no frequencies, so logomaker can't draw it
        empty = counts.index[counts.sum(axis=1) == 0]
        if len(empty):
            flank_len = len(counts)
            positions = [i - flank_len if side == "LEFT" else i + 1 for i in empty]
            logger.warning(
                f"Skipping {side.lower()} flank logo: position(s) {positions} "
                f"(relative to {mge_name}) are N in every sequence"
            )
            continue

        freqs = counts.div(counts.sum(axis=1), axis=0)

        fig, ax = plt.subplots(figsize=(10, 3))
        logomaker.Logo(freqs, ax=ax, shade_below=.5, fade_below=.5)
        ax.set_title(f"{side.capitalize()} Flank Sequence Logo")
        ax.set_xlabel("Position")
        ax.set_ylabel("Frequency")
        fig.tight_layout()
        fig.savefig(f"{output_dir}/graphs/{mge_name}_{side}_LOGO.png")
        plt.close(fig)
