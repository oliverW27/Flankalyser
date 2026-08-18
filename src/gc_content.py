"""
GC Content Analysis Module
Computes and visualizes GC content distributions for flanking sequences.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def gc_percent(seq):
    """Calculate GC percentage for a DNA sequence."""
    gc = seq.count('G') + seq.count('C')
    return 100 * gc / len(seq) if len(seq) > 0 else 0


def compute_gc_content(left_seqs, right_seqs, output_dir):
    """Compute GC content for unique left and right flanks and visualise distributions."""
    # Deduplicate sequence
    unique_left_seqs = list(set(left_seqs))
    unique_right_seqs = list(set(right_seqs))
    
    # Compute GC% for each unique sequence
    left_gc = [gc_percent(seq[:15]) for seq in unique_left_seqs]
    right_gc = [gc_percent(seq[15:]) for seq in unique_right_seqs]

    # Combine for total GC
    all_seqs = unique_left_seqs + unique_right_seqs
    all_gc_values = [gc_percent(seq) for seq in all_seqs]
    total_gc_mean = sum(all_gc_values) / len(all_gc_values) if all_gc_values else 0

    # Create DataFrame for plotting
    df = pd.DataFrame({
        'GC': left_gc + right_gc,
        'side': ['left'] * len(left_gc) + ['right'] * len(right_gc)
    })

    # Plot GC% distribution
    possible_gc_values = [i * 100 / 30 for i in range(31)]
    bins = possible_gc_values

    plt.figure(figsize=(8, 5))
    sns.histplot(
        data=df,
        x='GC',
        hue='side',
        bins=bins,
        palette=['#1f77b4', '#ff7f0e'],
        alpha=0.6
    )
    plt.title("GC Content Distribution of Flanks")
    plt.xlabel("GC %")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/graphs/gc_content_distribution.png")
    plt.close()

    # Compute summary statistics
    summary = df.groupby('side')['GC'].agg(['mean', 'median', 'min', 'max', 'count'])
    summary.loc['total'] = [total_gc_mean, None, None, None, len(all_gc_values)]

    # Save summary
    summary.to_csv(f"{output_dir}/files/gc_content_summary.tsv", sep='\t')
