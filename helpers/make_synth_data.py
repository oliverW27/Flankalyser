import random
from pathlib import Path

NUCLEOTIDES = ['A', 'T', 'G', 'C']

# Edit below to change the number of contexts or their length
N = 1000
FLANK_LEN = 50

# Random insertion - A,T,G,C at equal proportion at every position - no target site preference at all
def no_preference(input_file, flank_len, n):
    file = []
    for i in range(0,n):
        context = ''.join(random.choices(NUCLEOTIDES, k=flank_len*2))
        row = "\t".join(str(item) for item in [i, context[:flank_len], context[flank_len:]])
        file.append(row)
    with open(input_file, 'w') as f:
        for r in file:
            f.write(r + "\n")

# Preference to insert next to, or into a specific sequence at a specific proportion
def seq_preference(input_file, flank_len, n, specificity, seq,  side: str):
    if specificity > 100 or specificity < 0:
            raise ValueError("Must be an integer between 0 and 100, where 0 is only GC and 100 is only AT. 50 is ")
    file = []
    for i in range(0,n):
        # Element inserts downstream of motif - motif will be LEFT of the element
        if side == "l":
            context_l = ''.join(random.choices(NUCLEOTIDES, k=flank_len-len(seq)))
            context_r = ''.join(random.choices(NUCLEOTIDES, k=flank_len))
            if random.uniform(0,100) < specificity:
                context_l += seq
            else:
                context_l += ''.join(random.choices(NUCLEOTIDES, k=len(seq)))

        # Element inserts upstream of motif - motif will be RIGHT of the element
        elif side == "r":
            context_l = ''.join(random.choices(NUCLEOTIDES, k=flank_len))
            context_r = ''.join(random.choices(NUCLEOTIDES, k=flank_len-len(seq)))
            if random.uniform(0,100) < specificity:
                context_r = seq + context_r
            else:
                context_r += ''.join(random.choices(NUCLEOTIDES, k=len(seq)))

        # Element inserts into the motif - motif will be SURROUNDING the element.
        else:
            context = ''.join(random.choices(NUCLEOTIDES, k=flank_len*2 - len(seq)))
            if random.uniform(0,100) < specificity:
                context = context[0:flank_len-(len(seq)//2)] + seq + context[flank_len-(len(seq)//2):]
                context_l = context[0:flank_len]
                context_r = context[flank_len:]

            else:
                context += ''.join(random.choices(NUCLEOTIDES, k=len(seq)))

                context_l = context[0:flank_len]
                context_r = context[flank_len:]
        row = "\t".join(str(item) for item in [i, context_l, context_r])
        file.append(row)
    with open(input_file, 'w') as f:
        for r in file:
            f.write(r + "\n")

def nucl_preference(input_file, flank_len, n, specificity: int):
    if specificity > 100 or specificity < 0:
        raise ValueError("Must be an integer between 0 and 100, where 0 is only GC and 100 is only AT. 50 is ")
    file = []
    for i in range(0,n):
        context = "".join(random.choices(NUCLEOTIDES, weights=[specificity, specificity, 100-specificity, 100-specificity], k=flank_len*2))
        row = "\t".join(str(item) for item in [i, context[:flank_len], context[flank_len:]])
        file.append(row)
    with open(input_file, 'w') as f:
        for r in file:
            f.write(r + "\n")


if __name__ == '__main__':
    # Write to examples/toy_examples/ in the repo, wherever the script is run from
    file_path = Path(__file__).resolve().parent.parent / "examples" / "toy_examples"
    file_path.mkdir(parents=True, exist_ok=True)

    ### No preference at all
    no_preference(file_path / "flanks_no_preference.txt", FLANK_LEN, N)


    ###  Nucleotide preferences ###
    # Contexts are 20% AT
    nucl_preference(file_path / "flanks_low_AT.txt", FLANK_LEN, N, 20)
    # Contexts are 75% AT
    nucl_preference(file_path / "flanks_high_AT.txt", FLANK_LEN, N, 75)


    #### Sequence preferences ###
    # Inserts downstream of "CTAG" 80% of the time
    seq_preference(file_path / "flanks_strong_short_seq.txt", FLANK_LEN, N, 80, "CTAG", "l")
    # Inserts into "GCCGATCGCATCGGC" 10% of the time
    seq_preference(file_path / "flanks_weak_long_seq.txt", FLANK_LEN, N, 10, "GCCGATCGCATCGGC", "both")