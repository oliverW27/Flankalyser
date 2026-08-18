import matplotlib.pyplot as plt
from collections import Counter
import random

def detect_tsd_all_vs_all(
    flank_len,
    combined_props,
    left_seqs,
    right_seqs,
    output,
    min_len=2,
    max_len=15,
    compare_with_random=True,
    random_seed=42
):
    """Detect TSDs by comparing paired left/right sequences."""
    output_path = output + "/graphs/tsd_length_distribution.png"
    tsd_text_output = output + "/files/tsd_counts.txt"
    
    random.seed(random_seed)
    tsd_counter = Counter()
    
    # Compare paired sequences (same index)
    for i in range(min(len(left_seqs), len(right_seqs))):
        left = left_seqs[i]
        right = right_seqs[i]
        
        # Check for TSD: end of left matches start of right
        max_check = min(len(left), len(right), max_len)
        for k in range(max_check, min_len - 1, -1):
            if left[-k:] == right[:k]:
                tsd_counter[left[-k:]] += 1
                break  # Take longest match only
    
    # Write TSD sequences and counts
    with open(tsd_text_output, "w") as f:
        f.write("TSD_Sequence\tCount\tLength\n")
        for tsd, count in tsd_counter.most_common():
            f.write(f"{tsd}\t{count}\t{len(tsd)}\n")
    
    # Get length distribution
    observed_lengths = [len(tsd) for tsd in tsd_counter.elements()]
    obs_len_counts = Counter(observed_lengths)
    
    # Random comparison: shuffle the sequences
    if compare_with_random:
        rand_counter = Counter()
        
        # Shuffle right sequences to break real TSDs
        shuffled_right = right_seqs.copy()
        random.shuffle(shuffled_right)
        
        for i in range(min(len(left_seqs), len(shuffled_right))):
            left = left_seqs[i]
            right = shuffled_right[i]
            
            max_check = min(len(left), len(right), max_len)
            for k in range(max_check, min_len - 1, -1):
                if left[-k:] == right[:k]:
                    rand_counter[left[-k:]] += 1
                    break
        
        rand_lengths = [len(tsd) for tsd in rand_counter.elements()]
        rand_len_counts = Counter(rand_lengths)
    else:
        rand_len_counts = {}
    
    # Plot length distribution
    all_lengths = sorted(set(obs_len_counts.keys()) | set(rand_len_counts.keys()))
    obs_vals = [obs_len_counts.get(l, 0) for l in all_lengths]
    rand_vals = [rand_len_counts.get(l, 0) for l in all_lengths]
    
    plt.figure(figsize=(10, 6))
    if compare_with_random:
        plt.bar([l - 0.2 for l in all_lengths], obs_vals, width=0.4, 
                color='tab:blue', label='Observed', alpha=0.7)
        plt.bar([l + 0.2 for l in all_lengths], rand_vals, width=0.4,
                color='tab:orange', label='Shuffled', alpha=0.7)
    else:
        plt.bar(all_lengths, obs_vals, color='tab:blue', label='Observed', alpha=0.7)
    
    plt.xlabel("TSD Length", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.title("Target Site Duplication Length Distribution", fontsize=16)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    
    print(f"Found {len(tsd_counter)} unique TSDs")
    print(f"Total TSD events: {sum(tsd_counter.values())}")
    
    return len(tsd_counter)