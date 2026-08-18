"""
K-mer Count Analysis Module
Analyses unique k-mer distributions in flanking sequences vs random controls.
"""

import matplotlib.pyplot as plt
import logging
from src.utils import generate_position_weighted_random_seq

logger = logging.getLogger(__name__)


def get_unique_kmers(seqs, k):
    """Return a set of all unique k-mers from a list of sequences."""
    kmers = set()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            kmers.add(kmer)
    logger.debug(f"Extracted {len(kmers)} unique {k}-mers from {len(seqs)} sequences")
    return kmers


def plot_unique_kmer_counts_both(left_seqs, right_seqs, combined_props, flank_len,
                                min_len=4, n_random=1):
    """Plot number of unique k-mers from length min_len to max_len in actual vs. random sequences,
    for both left and right flanks side-by-side.

    Parameters:
    - left_seqs, right_seqs: lists of sequences
    - combined_props: dict with keys 'A','C','G','T' and 60-length lists (first 30 left, last 30 right)
    - flank_len: length of flank region
    - min_len: minimum k-mer length
    - n_random: number of random sequences generated per original sequence
    """
    ks = list(range(min_len, flank_len + 1))
    output_path = "graphs/kmer_distribution.png"

    logger.info(
        f"Starting unique k-mer count analysis "
        f"(k range: {min_len}–{flank_len}, n_random={n_random})"
    )

    def calc_counts(seqs, flank):
        """Calculate actual vs expected k-mer counts for given sequences."""
        actual_counts = []
        expected_counts = []
        
        logger.info(f"Processing {len(seqs)} sequences for {flank} flank...")

        for k in ks:
            logger.debug(f"Counting {k}-mers for {flank} flank")

            # Count actual unique k-mers
            actual_kmers = get_unique_kmers(seqs, k)
            actual_counts.append(len(actual_kmers))

            # Generate random k-mers using position-weighted probabilities
            random_kmers = set()
            for seq in seqs:
                for _ in range(n_random):
                    if flank == "left":
                        random_seq = generate_position_weighted_random_seq(combined_props, 0, flank_len)
                    elif flank == "right":
                        random_seq = generate_position_weighted_random_seq(combined_props, flank_len + 1, 2*flank_len + 1)
                    random_kmers.update(get_unique_kmers([random_seq], k))

            expected_counts.append(len(random_kmers))

            logger.debug(
                f"{flank.capitalize()} flank {k}-mers: "
                f"{len(actual_kmers)} actual, {len(random_kmers)} expected"
            )
            
        return actual_counts, expected_counts

    # Calculate counts for both flanks
    left_actual, left_expected = calc_counts(left_seqs, 'left')
    right_actual, right_expected = calc_counts(right_seqs, 'right')

    logger.info("Finished calculating unique k-mer counts, generating plots...")

    # Create side-by-side plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    # Left flank plot
    axes[0].plot(ks, left_actual, label='Actual', marker='o', color='tab:blue')
    axes[0].plot(ks, left_expected, label='Expected (random)', marker='s', color='tab:orange')
    axes[0].set_title("Left Flank Unique k-mer Counts")
    axes[0].set_xlabel("k-mer Length")
    axes[0].set_ylabel("Unique k-mers")
    axes[0].grid(True)
    axes[0].legend()

    # Right flank plot
    axes[1].plot(ks, right_actual, label='Actual', marker='o', color='tab:blue')
    axes[1].plot(ks, right_expected, label='Expected (random)', marker='s', color='tab:orange')
    axes[1].set_title("Right Flank Unique k-mer Counts")
    axes[1].set_xlabel("k-mer Length")
    axes[1].grid(True)
    axes[1].legend()

    plt.tight_layout()

    try:
        plt.savefig(output_path)
        logger.info(f"Saved k-mer distribution plot to {output_path}")
    except Exception as e:
        logger.exception(f"Failed to save k-mer plot: {e}")
    finally:
        plt.close()
