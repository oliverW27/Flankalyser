import random
import numpy as np
import matplotlib.pyplot as plt
from src.utils import generate_position_weighted_random_seq

def reverse_complement(left, right):
    """Check if left is reverse complement of right."""
    return left.translate(str.maketrans("ACGT", "TGCA"))[::-1] == right

def generate_random_seq(length=0):
    return ''.join(random.choices("ACGT", k=length))

def plot_palindrome_results(results, random_results, output_path="graphs/palindrome_centered_barplot.png"):
    ks = sorted(results.keys())
    actual_counts = [results[k] for k in ks]
    random_counts = [random_results[k] for k in ks]

    x = np.arange(len(ks))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8,6))
    ax.bar(x - width/2, actual_counts, width, label='Actual', color='tab:blue')
    ax.bar(x + width/2, random_counts, width, label='Random', color='tab:orange')

    ax.set_xticks(x)
    ax.set_xticklabels(ks)
    ax.set_xlabel("Palindrome Length")
    ax.set_ylabel("Count of Palindromes")
    ax.set_title("Longest Centered Palindromes Across Junction")
    ax.legend()
    ax.grid(True, axis='y')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
def process_summary_file_centered_palindromes(left_seqs, right_seqs, combined_props, path, min_len=5, max_len=8):
    """
    Corrected function to process summary file for centered palindromes.
    The left_seqs and right_seqs arguments are accepted to match the original
    function call but are no longer used in the logic.
    """
    results = {k: 0 for k in range(min_len, max_len + 1)}
    random_results = {k: 0 for k in range(min_len, max_len + 1)}

    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            
            accession, lefts_raw, rights_raw = parts[:3]

            if lefts_raw == "none" or rights_raw == "none":
                continue

            lefts = lefts_raw.split("|")
            rights = rights_raw.split("|")

            # Assume all flank lengths are consistent based on the first sequence
            flank_len = len(lefts[0])

            # Iterate through actual pairs from the file
            for left in lefts:
                for right in rights:
                    # Check the actual pair
                    for k in range(min_len, max_len + 1):
                        if len(left) >= k and len(right) >= k and reverse_complement(left[-k:], right[:k]):
                            results[k] += 1
                    
                    # For each actual comparison, perform one random comparison
                    rand_left = generate_position_weighted_random_seq(combined_props, 0, flank_len)
                    rand_right = generate_position_weighted_random_seq(combined_props, flank_len + 1, (flank_len * 2) + 1)
                    
                    for k in range(min_len, max_len + 1):
                        if reverse_complement(rand_left[-k:], rand_right[:k]):
                            random_results[k] += 1

    plot_palindrome_results(results, random_results)