import random
from collections import defaultdict


def dinucleotide_shuffle(seq, max_attempts=200):
    """Shuffle a sequence while preserving its dinucleotide frequencies (Altschul & Erickson 1985).

    Defined once here so every module that needs a composition-preserving
    shuffled-sequence null (palindrome detection, melting profile controls,
    etc.) calls the same implementation instead of each re-implementing it.
    """
    seq = seq.upper()
    for _ in range(max_attempts):
        edges = defaultdict(list)
        for i in range(len(seq) - 1):
            edges[seq[i]].append(seq[i + 1])
        for base in edges:
            random.shuffle(edges[base])
        result = [seq[0]]
        current = seq[0]
        while edges[current]:
            next_base = edges[current].pop()
            result.append(next_base)
            current = next_base
        if len(result) == len(seq):
            return "".join(result)
    raise RuntimeError(f"Dinucleotide shuffle failed after {max_attempts} attempts")


def generate_position_weighted_random_seq(prop_dict, start, end):
    """
    Generate a random DNA sequence using position-specific nucleotide probabilities.

    Args:
        prop_dict (dict): Dictionary with keys 'A', 'C', 'G', 'T' mapping to lists of probabilities.
                          Each list must have length >= end.
        start (int): Start index (inclusive) in the position list.
        end (int): End index (exclusive) in the position list.

    Returns:
        str: Random DNA sequence of length (end - start) where each base is chosen based on the
             nucleotide proportions at that position.
    """
    sequence = []
    for i in range(start, end):
        try:
            weights = [prop_dict[nt][i] for nt in "ACGT"]
        except (KeyError, IndexError) as e:
            raise ValueError(f"Invalid prop_dict structure or position index out of range at position {i}.") from e
        base = random.choices("ACGT", weights=weights, k=1)[0]
        sequence.append(base)
    return ''.join(sequence)