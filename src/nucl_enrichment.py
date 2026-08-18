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


def _count_trinucleotides(seqs):
    """Count all trinucleotides across a list of sequences."""
    counts = Counter()
    total = 0
    for seq in seqs:
        for i in range(len(seq) - 2):
            counts[seq[i:i+3]] += 1
            total += 1
    return counts, total


def analyse_dinucleotide_enrichment_weighted(left_seqs, right_seqs, baseline_props, output_dir, flank_len, inner_bp=7):
    """
    Dinucleotide enrichment: inner N bp vs outer (flank_len - N) bp of each flank.
    Inner = closer to insertion site, outer = further away.

    NOTE: baseline_props is unused now but kept for interface compatibility.
    """
    split = flank_len - inner_bp

    for flank, seqs in [('left', left_seqs), ('right', right_seqs)]:
        if flank == 'left':
            # Inner = last 25bp (near insertion), outer = first 25bp
            inner_seqs = [s[split:] for s in seqs]
            outer_seqs = [s[:split] for s in seqs]
        else:
            # Inner = first 25bp (near insertion), outer = last 25bp
            inner_seqs = [s[:split] for s in seqs]
            outer_seqs = [s[split:] for s in seqs]

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

        df.to_csv(f"{output_dir}/files/{flank}_dinucleotide_enrichment.tsv", sep='\t')

        plt.figure(figsize=(6, 5))
        sns.heatmap(df.astype(float), cmap='coolwarm', center=0, annot=True, fmt=".2f")
        plt.title(f"{flank.capitalize()} Flank Dinucleotide Enrichment\n(log2 inner/outer)")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/graphs/{flank}_dinucleotide_enrichment.png")
        plt.close()


def analyse_trinucleotide_enrichment_weighted(left_seqs, right_seqs, baseline_props, output_dir, flank_len, inner_bp=25, top_n=20):
    """
    Trinucleotide enrichment: inner N bp vs outer (flank_len - N) bp of each flank.

    NOTE: baseline_props is unused now but kept for interface compatibility.
    """
    split = flank_len - inner_bp

    for flank, seqs in [('left', left_seqs), ('right', right_seqs)]:
        if flank == 'left':
            inner_seqs = [s[split:] for s in seqs]
            outer_seqs = [s[:split] for s in seqs]
        else:
            inner_seqs = [s[:split] for s in seqs]
            outer_seqs = [s[split:] for s in seqs]

        inner_counts, inner_total = _count_trinucleotides(inner_seqs)
        outer_counts, outer_total = _count_trinucleotides(outer_seqs)

        enrichment = {}
        for a in 'ACGT':
            for b in 'ACGT':
                for c in 'ACGT':
                    trinuc = a + b + c
                    inner_freq = (inner_counts.get(trinuc, 0) + 0.5) / (inner_total + 32)
                    outer_freq = (outer_counts.get(trinuc, 0) + 0.5) / (outer_total + 32)
                    enrichment[trinuc] = np.log2(inner_freq / outer_freq)

        sorted_enrichment = sorted(enrichment.items(), key=lambda x: -abs(x[1]))[:top_n]
        df = pd.DataFrame(sorted_enrichment, columns=['motif', 'log2_enrichment'])

        df.to_csv(f"{output_dir}/files/{flank}_trinucleotide_enrichment.tsv", sep='\t', index=False)

        plt.figure(figsize=(10, 4))
        sns.barplot(data=df, x='motif', y='log2_enrichment', palette='viridis')
        plt.axhline(0, color='black', linestyle='--')
        plt.xticks(rotation=45)
        plt.title(f"{flank.capitalize()} Trinucleotide Enrichment\n(log2 inner/outer)")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/graphs/{flank}_trinucleotide_enrichment.png")
        plt.close()