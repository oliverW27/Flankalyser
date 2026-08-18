"""
Sequence Logo Generation Module
Creates sequence logos for flanking sequences using logomaker.
"""
import pandas as pd
import matplotlib.pyplot as plt
import logomaker


def generate_sequence_logo(left_seqs, right_seqs, output_dir):
    """Generate sequence logos for left and right flanks.
    
    Parameters:
    - left_seqs: list of left flank sequences (assumed equal length)
    - right_seqs: list of right flank sequences (assumed equal length)
    - output_dir: directory to save the logo images
    """
    
    def seqs_to_counts(seqs):
        """Convert list of sequences to position frequency matrix (DataFrame)."""
        length = len(seqs[0])
        counts = {base: [0] * length for base in 'ACGT'}
        
        for seq in seqs:
            for i, nt in enumerate(seq):
                if nt in counts:
                    counts[nt][i] += 1
                    
        return pd.DataFrame(counts).T  # Rows = bases, columns = positions
    
    # Filter sequences of the same length (in case input is variable)
    left_seqs = [s for s in left_seqs if len(s) == len(left_seqs[0])]
    right_seqs = [s for s in right_seqs if len(s) == len(right_seqs[0])]
    
    # Build count matrices
    left_counts = seqs_to_counts(left_seqs).T  # now columns = A,C,G,T; rows = positions
    right_counts = seqs_to_counts(right_seqs).T
    
    # Normalize to frequencies
    left_freq = left_counts.div(left_counts.sum(axis=1), axis=0)
    right_freq = right_counts.div(right_counts.sum(axis=1), axis=0)
    
    # Plot left flank logo
    plt.figure(figsize=(10, 3))
    logo_left = logomaker.Logo(left_freq, shade_below=.5, fade_below=.5)
    plt.title("Left Flank Sequence Logo")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/graphs/left_flank_sequence_logo.png")
    plt.close()
    
    # Plot right flank logo
    plt.figure(figsize=(10, 3))
    logo_right = logomaker.Logo(right_freq, shade_below=.5, fade_below=.5)
    plt.title("Right Flank Sequence Logo")
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/graphs/right_flank_sequence_logo.png")
    plt.close()