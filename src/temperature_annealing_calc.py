#!/usr/bin/env python3
"""
Melting analysis of IS element flanking sequences.

Three per-position profiles centered on the insertion site:
  1. Mean hydrogen bond count (2 = AT, 3 = GC)
  2. Nearest-neighbor ΔG (SantaLucia 1998 unified params)
  3. Peyrard-Bishop model breathing probability

Usage: python melting_analysis.py input.tsv [output_prefix] [max_pb_seqs]

Input: tab-separated (sample_id  left_flank  right_flank)
  Oriented as: LEFT_FLANK | IS | RIGHT_FLANK
  Ns are skipped; all-N sequences excluded; flanks deduplicated.

References:
  SantaLucia J Jr (1998) PNAS 95:1460-1465  (NN params)
  Peyrard M & Bishop AR (1989) Phys Rev Lett 62:2755-2758  (PB model)
  Dauxois T, Peyrard M, Bishop AR (1993) Phys Rev E 47:R44  (anharmonic stacking)
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

# --- SantaLucia (1998) unified NN ΔG°37 in kcal/mol ---
NN_DG = {
    'AA': -1.00, 'TT': -1.00, 'AT': -0.88, 'TA': -0.58,
    'CA': -1.45, 'TG': -1.45, 'GT': -1.44, 'AC': -1.44,
    'CT': -1.28, 'AG': -1.28, 'GA': -1.30, 'TC': -1.30,
    'CG': -2.17, 'GC': -2.24, 'GG': -1.84, 'CC': -1.84,
}
HBONDS = {'A': 2, 'T': 2, 'G': 3, 'C': 3}


def parse_input(filepath):
    """Read TSV, split | flanks, deduplicate, skip all-N."""
    lefts, rights = set(), set()
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            # Split pipe-separated flanks into independent sequences
            for lf in parts[1].upper().split('|'):
                lf = lf.strip()
                if lf and set(lf) != {'N'}:
                    lefts.add(lf)
            for rf in parts[2].upper().split('|'):
                rf = rf.strip()
                if rf and set(rf) != {'N'}:
                    rights.add(rf)
    return list(lefts), list(rights)


def profile_with_sem(seqs, per_seq_fn):
    """Mean and SEM across per-sequence profiles."""
    profiles = [per_seq_fn(s) for s in seqs if per_seq_fn(s) is not None]
    if not profiles:
        return None, None
    max_len = max(len(p) for p in profiles)
    arr = np.full((len(profiles), max_len), np.nan)
    for i, p in enumerate(profiles):
        arr[i, :len(p)] = p
    mean = np.nanmean(arr, axis=0)
    sem = np.nanstd(arr, axis=0) / np.sqrt(np.sum(~np.isnan(arr), axis=0))
    return mean, sem


def hbond_per_seq(seq):
    return np.array([HBONDS.get(b, np.nan) for b in seq])


def dg_per_seq(seq):
    vals = [NN_DG.get(seq[i:i+2], np.nan) for i in range(len(seq) - 1)]
    return np.array(vals)


def pb_breathing_profile(seqs, T=310, max_seqs=None, pad=10):
    """
    Peyrard-Bishop-Dauxois breathing probability per position.
    Uses transfer-integral forward-backward with G-padding to suppress
    edge artefacts. Returns (mean, sem) arrays.
    """
    kB = 8.617e-5   # eV/K
    beta = 1.0 / (kB * T)

    # Morse potential depth (eV) — AT weaker (2 H-bonds), GC stronger (3)
    D_map = {'A': 0.05, 'T': 0.05, 'G': 0.075, 'C': 0.075}
    alpha = 4.2      # Morse width (Å⁻¹)
    K_s = 0.025      # stacking spring constant (eV/Å²)
    rho = 0.5        # stacking anharmonicity amplitude
    delta_s = 0.35   # stacking anharmonicity range (Å⁻¹)
    threshold = 0.5  # opening threshold (Å)

    n_grid = 200
    y = np.linspace(-0.3, 5.0, n_grid)
    dy = y[1] - y[0]
    open_mask = y > threshold

    # Precompute symmetric stacking kernel (sequence-independent in basic PBD)
    Y1, Y2 = np.meshgrid(y, y, indexing='ij')
    W = 0.5 * K_s * (1.0 + rho * np.exp(-delta_s * (Y1 + Y2))) * (Y1 - Y2) ** 2
    Kmat = np.exp(-beta * W)

    # Filter: no Ns, consistent length
    valid = [s for s in seqs if 'N' not in s]
    if not valid:
        return None, None
    lengths = [len(s) for s in valid]
    target_len = max(set(lengths), key=lengths.count)
    valid = [s for s in valid if len(s) == target_len]
    if not valid:
        return None, None

    if max_seqs and len(valid) > max_seqs:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(valid), max_seqs, replace=False)
        valid = [valid[i] for i in idx]

    pad_seq = 'G' * pad  # maximally stable padding

    all_probs = np.zeros((len(valid), target_len))
    for si, seq in enumerate(valid):
        padded = pad_seq + seq + pad_seq
        N = len(padded)

        # Morse potentials per position
        V = np.zeros((N, n_grid))
        for i, base in enumerate(padded):
            d = D_map.get(base, 0.0625)
            V[i] = d * (np.exp(-alpha * y) - 1.0) ** 2

        # Forward pass: f_n(y) = exp(-β V_n) * ∫ K(y',y) f_{n-1}(y') dy'
        fwd = np.zeros((N, n_grid))
        fwd[0] = np.exp(-beta * V[0])
        fwd[0] /= (np.sum(fwd[0]) * dy)
        for n in range(1, N):
            conv = Kmat @ fwd[n - 1] * dy
            fwd[n] = np.exp(-beta * V[n]) * conv
            s = np.sum(fwd[n]) * dy
            if s > 0:
                fwd[n] /= s

        # Backward pass: b_n(y) = ∫ K(y,y') exp(-β V_{n+1}(y')) b_{n+1}(y') dy'
        bwd = np.zeros((N, n_grid))
        bwd[N - 1] = np.ones(n_grid)
        bwd[N - 1] /= (np.sum(bwd[N - 1]) * dy)
        for n in range(N - 2, -1, -1):
            temp = np.exp(-beta * V[n + 1]) * bwd[n + 1]
            conv = Kmat @ temp * dy
            s = np.sum(conv) * dy
            if s > 0:
                conv /= s
            bwd[n] = conv

        # Marginal opening probabilities, trimming padding
        for n in range(pad, pad + target_len):
            marginal = fwd[n] * bwd[n]
            total = np.sum(marginal) * dy
            if total > 0:
                all_probs[si, n - pad] = np.sum(marginal[open_mask]) * dy / total

        if (si + 1) % 500 == 0:
            print(f"  PB: {si+1}/{len(valid)}", file=sys.stderr)

    mean = np.mean(all_probs, axis=0)
    sem = np.std(all_probs, axis=0) / np.sqrt(len(valid))
    return mean, sem


def shuffle_seqs(seqs, rng):
    """Shuffle each sequence to destroy positional signal, keep composition."""
    out = []
    for s in seqs:
        arr = list(s)
        rng.shuffle(arr)
        out.append(''.join(arr))
    return out


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    infile = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "melting_analysis"
    max_pb = int(sys.argv[3]) if len(sys.argv) > 3 else 5000

    print("Parsing...", file=sys.stderr)
    lefts, rights = parse_input(infile)
    print(f"  {len(lefts)} unique left, {len(rights)} unique right flanks", file=sys.stderr)

    rng = np.random.default_rng(42)
    shuf_lefts = shuffle_seqs(lefts, rng)
    shuf_rights = shuffle_seqs(rights, rng)

    left_len = len(lefts[0]) if lefts else 0
    right_len = len(rights[0]) if rights else 0

    # --- Compute all profiles ---
    print("H-bond profiles...", file=sys.stderr)
    hb_l, hb_l_sem = profile_with_sem(lefts, hbond_per_seq)
    hb_r, hb_r_sem = profile_with_sem(rights, hbond_per_seq)
    hb_sl, _ = profile_with_sem(shuf_lefts, hbond_per_seq)
    hb_sr, _ = profile_with_sem(shuf_rights, hbond_per_seq)

    print("NN ΔG profiles...", file=sys.stderr)
    dg_l, dg_l_sem = profile_with_sem(lefts, dg_per_seq)
    dg_r, dg_r_sem = profile_with_sem(rights, dg_per_seq)
    dg_sl, _ = profile_with_sem(shuf_lefts, dg_per_seq)
    dg_sr, _ = profile_with_sem(shuf_rights, dg_per_seq)

    print(f"PB breathing (max {max_pb} seqs/side)...", file=sys.stderr)
    pb_l, pb_l_sem = pb_breathing_profile(lefts, max_seqs=max_pb)
    pb_r, pb_r_sem = pb_breathing_profile(rights, max_seqs=max_pb)
    pb_sl, _ = pb_breathing_profile(shuf_lefts, max_seqs=max_pb)
    pb_sr, _ = pb_breathing_profile(shuf_rights, max_seqs=max_pb)

    # --- X-axes relative to insertion ---
    x_l = np.arange(-left_len, 0)
    x_r = np.arange(1, right_len + 1)
    x_l_dg = x_l[:-1] + 0.5   # dinucleotide step centers
    x_r_dg = x_r[:-1] + 0.5

    # --- Plot ---
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    dk = dict(linewidth=1.0)
    sk = dict(linewidth=1.0, linestyle='--', alpha=0.4, color='grey')
    fa = 0.12

    # Panel 1: H-bonds
    ax = axes[0]
    ax.plot(x_l, hb_l, color='#2166ac', **dk, label='Left flank')
    ax.fill_between(x_l, hb_l - hb_l_sem, hb_l + hb_l_sem, color='#2166ac', alpha=fa)
    ax.plot(x_r, hb_r, color='#b2182b', **dk, label='Right flank')
    ax.fill_between(x_r, hb_r - hb_r_sem, hb_r + hb_r_sem, color='#b2182b', alpha=fa)
    if hb_sl is not None:
        ax.plot(x_l, hb_sl, **sk, label='Shuffled')
        ax.plot(x_r, hb_sr, **sk)
    ax.axvline(0.5, color='grey', ls='--', alpha=0.4)
    ax.set_ylabel('Mean H-bonds per bp')
    ax.set_title('Hydrogen bond count (2 = AT, 3 = GC)')
    ax.legend(fontsize=8, loc='upper right')

    # Panel 2: NN ΔG  (less negative = easier to melt = higher on plot)
    ax = axes[1]
    ax.plot(x_l_dg, dg_l, color='#2166ac', **dk)
    ax.fill_between(x_l_dg, dg_l - dg_l_sem, dg_l + dg_l_sem, color='#2166ac', alpha=fa)
    ax.plot(x_r_dg, dg_r, color='#b2182b', **dk)
    ax.fill_between(x_r_dg, dg_r - dg_r_sem, dg_r + dg_r_sem, color='#b2182b', alpha=fa)
    if dg_sl is not None:
        ax.plot(x_l_dg, dg_sl, **sk)
        ax.plot(x_r_dg, dg_sr, **sk)
    ax.axvline(0.5, color='grey', ls='--', alpha=0.4)
    ax.set_ylabel('ΔG°₃₇ (kcal/mol)')
    ax.set_title('Nearest-neighbor stacking ΔG (SantaLucia 1998)')
    ax.annotate('↑ easier to melt', xy=(0.01, 0.95), xycoords='axes fraction',
                fontsize=7, color='grey', va='top')

    # Panel 3: PB breathing
    ax = axes[2]
    if pb_l is not None:
        ax.plot(x_l, pb_l, color='#2166ac', **dk)
        ax.fill_between(x_l, pb_l - pb_l_sem, pb_l + pb_l_sem, color='#2166ac', alpha=fa)
    if pb_r is not None:
        ax.plot(x_r, pb_r, color='#b2182b', **dk)
        ax.fill_between(x_r, pb_r - pb_r_sem, pb_r + pb_r_sem, color='#b2182b', alpha=fa)
    if pb_sl is not None:
        ax.plot(x_l, pb_sl, **sk)
    if pb_sr is not None:
        ax.plot(x_r, pb_sr, **sk)
    ax.axvline(0.5, color='grey', ls='--', alpha=0.4)
    ax.set_ylabel('P(opening)')
    ax.set_title('Peyrard-Bishop breathing probability (T = 310 K)')
    ax.set_xlabel('Position relative to IS insertion site')

    plt.tight_layout()
    plt.savefig(f"{prefix}.png", dpi=200, bbox_inches='tight')
    plt.savefig(f"{prefix}.pdf", bbox_inches='tight')
    print(f"Saved {prefix}.png / .pdf", file=sys.stderr)


if __name__ == '__main__':
    main()
