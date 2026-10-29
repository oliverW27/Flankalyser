import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


def _count_dinucleotides(seqs):
    """Count all dinucleotides across a list of sequences."""
    counts = Counter()
    total = 0
    for seq in seqs:
        for i in range(len(seq) - 1):
            counts[seq[i:i+2]] += 1
            total += 1
    return counts, total


def analyse_dinucleotide_enrichment_weighted(left_seqs, right_seqs, output_dir, mge_name, inner_bp=10):
    """
    Dinucleotide enrichment: the inner_bp nearest the insertion site vs. the
    inner_bp farthest from it -- symmetric, matched-size windows, so the
    middle of the flank (between the two windows) is left out of both.
    Inner = closer to insertion site, outer = farthest away.
    """
    dfs = {}
    for flank, seqs in [('left', left_seqs), ('right', right_seqs)]:
        if flank == 'left':
            # Inner = last inner_bp (near insertion), outer = first inner_bp (far)
            inner_seqs = [s[-inner_bp:] for s in seqs]
            outer_seqs = [s[:inner_bp] for s in seqs]
        else:
            # Inner = first inner_bp (near insertion), outer = last inner_bp (far)
            inner_seqs = [s[:inner_bp] for s in seqs]
            outer_seqs = [s[-inner_bp:] for s in seqs]

        inner_counts, inner_total = _count_dinucleotides(inner_seqs)
        outer_counts, outer_total = _count_dinucleotides(outer_seqs)

        # Log2 enrichment of inner freq vs outer freq
        df = pd.DataFrame(index=list('ACGT'), columns=list('ACGT'), dtype=float)
        for a in 'ACGT':
            for b in 'ACGT':
                dinuc = a + b
                inner_freq = (inner_counts.get(dinuc, 0) + 0.5) / (inner_total + 8)  # pseudocount
                outer_freq = (outer_counts.get(dinuc, 0) + 0.5) / (outer_total + 8)
                df.loc[a, b] = np.log2(inner_freq / outer_freq)

        df.to_csv(f"{output_dir}/files/{mge_name}_{flank.upper()}_DINUC_ENRICHMENT.tsv", sep='\t')
        dfs[flank] = df

    # One figure, left flank and right flank side by side, sharing a colour scale.
    vmax = max(dfs['left'].abs().values.max(), dfs['right'].abs().values.max())
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, flank in zip(axes, ('left', 'right')):
        sns.heatmap(dfs[flank], cmap='coolwarm', center=0, vmin=-vmax, vmax=vmax,
                    annot=True, fmt=".2f", ax=ax, cbar=(flank == 'right'))
        ax.set_title(f"{flank.capitalize()} flank")
    fig.suptitle(f"{mge_name} dinucleotide enrichment "
                 f"(log2 inner/outer, {inner_bp}bp nearest vs. {inner_bp}bp farthest)")
    fig.tight_layout()
    fig.savefig(f"{output_dir}/graphs/{mge_name}_DINUC_ENRICH.png")
    plt.close(fig)
