import random

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