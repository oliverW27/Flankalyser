"""
Palindrome / inverted-repeat detection module.

Standalone script: reads a trimmed flanks TSV (acc_name\\tleft\\tright),
concatenates left+right per accession, searches for inverted repeats with
optional spacer and mismatches within a configurable window, compares
against a dinucleotide-shuffled null, and writes a summary + two plots:
  - null distribution (observed hits vs. shuffle distribution, with p-value)
  - arm length distribution (observed vs. null, averaged per shuffle)

Usage:
    python palindrome_analysis.py --input flanks.tsv --mge_name IS10L --outdir results
"""

import argparse
import csv
import logging
import random
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.utils import dinucleotide_shuffle

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


# ---------------------------------------------------------------------------
# Inverted repeat search
# ---------------------------------------------------------------------------

def _is_complement(a, b):
    return _COMPLEMENT.get(a) == b


def find_longest_inverted_repeat(seq, max_spacer, min_arm=3, max_mismatches=0,
                                 search_start=None, search_end=None,
                                 fixed_centre=None):
    """Find the longest inverted repeat (arms are reverse complements) with a spacer of 0+ nt.

    Optionally restrict the search to a window [search_start, search_end).
    If fixed_centre is set, only that index is tried as the palindrome centre.
    Returns dict with arm_len, left_start, right_end, spacer_size, arm_seq, mismatches.
    """
    seq = seq.upper()
    n = len(seq)
    if search_start is None:
        search_start = 0
    if search_end is None:
        search_end = n
    search_start = max(0, search_start)
    search_end = min(n, search_end)

    best = {"arm_len": 0, "left_start": None, "right_end": None,
            "spacer_size": None, "arm_seq": None, "mismatches": 0}

    # When fixed_centre is set, only test that one position.
    centres = [fixed_centre] if fixed_centre is not None else range(search_start, search_end)
    for centre in centres:
        for spacer_size in range(0, max_spacer + 1):
            left_end = centre - spacer_size // 2
            right_start = centre + (spacer_size + 1) // 2
            if left_end < search_start or right_start >= search_end:
                break
            arm_len = 0
            mismatches = 0
            l, r = left_end - 1, right_start
            while l >= search_start and r < search_end:
                if _is_complement(seq[l], seq[r]):
                    arm_len += 1
                elif mismatches < max_mismatches:
                    mismatches += 1
                    arm_len += 1
                else:
                    break
                l -= 1
                r += 1
            if arm_len >= min_arm and arm_len > best["arm_len"]:
                best = {
                    "arm_len": arm_len,
                    "left_start": l + 1,
                    "right_end": r,
                    "spacer_size": spacer_size,
                    "arm_seq": seq[l + 1:r],
                    "mismatches": mismatches,
                }
    return best


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _sig_stars(p):
    """Conventional significance stars for a one-sided p-value."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_null_distribution(null_hits, observed_hits, p_val, mge_name, output_path, enrichment=None):
    """Histogram of per-shuffle hit counts, with a line marking the observed value.

    This is the plot that makes the p-value visual: if the observed line sits
    off to the right of the null histogram, the enrichment is real.
    """
    stars = _sig_stars(p_val)
    fig, ax = plt.subplots(figsize=(8, 5))
    lo, hi = min(null_hits + [observed_hits]), max(null_hits + [observed_hits])
    bins = range(lo, hi + 2)
    ax.hist(null_hits, bins=bins, color="tab:orange", edgecolor="white",
             align="left", label="Null (shuffled)")
    ax.axvline(observed_hits, color="tab:blue", linewidth=2.5,
               label=f"Observed ({observed_hits})")
    ax.set_xlabel("Sequences hitting the pattern")
    ax.set_ylabel("Number of shuffles")
    ax.set_title(f"{mge_name}: observed vs. null hit distribution (p = {p_val:.4f} {stars})")
    annotation = f"p = {p_val:.4f} {stars}"
    if enrichment is not None:
        annotation = f"Enrichment: {enrichment:.2f}x\n" + annotation
    ax.text(0.98, 0.95, annotation, transform=ax.transAxes, ha="right", va="top",
            fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_arm_length_distribution(obs_arm_lens, null_arm_lens, n_shuffles, mge_name, output_path,
                                  p_val=None, enrichment=None):
    """Bar chart of arm-length distribution, observed vs. null.

    Null counts are averaged per shuffle so they're on the same scale as the
    single observed dataset (raw null counts are summed across all shuffles).

    p_val/enrichment (if given) are the overall hit-rate significance from
    plot_null_distribution's test, annotated here for reference -- not a
    separate per-arm-length test.
    """
    all_arms = sorted(set(obs_arm_lens) | set(null_arm_lens))
    if not all_arms:
        return
    obs_counts = [obs_arm_lens.count(a) for a in all_arms]
    null_counts = [null_arm_lens.count(a) / n_shuffles for a in all_arms]

    x = range(len(all_arms))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([i - width / 2 for i in x], obs_counts, width, label="Observed", color="tab:blue")
    ax.bar([i + width / 2 for i in x], null_counts, width,
           label="Null (avg per shuffle)", color="tab:orange")
    ax.set_xticks(list(x))
    ax.set_xticklabels(all_arms)
    ax.set_xlabel("Arm length (bp)")
    ax.set_ylabel("Number of hits")
    ax.set_title(f"{mge_name}: inverted repeat arm length, observed vs. null")
    if p_val is not None:
        stars = _sig_stars(p_val)
        annotation = f"p = {p_val:.4f} {stars}"
        if enrichment is not None:
            annotation = f"Arm-len enrichment: {enrichment:.2f}x\n" + annotation
        ax.text(0.98, 0.95, annotation, transform=ax.transAxes, ha="right", va="top",
                fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ---------------------------------------------------------------------------
# Input loading
def load_sequences(path, min_len=20):
    """Load (acc, left, right) from TSV; strip trailing Ns from right flank, dedupe.

    Returns list of (acc, concat_seq).
    """
    seqs = []
    seen = set()
    with open(path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue
            acc, left, right = row[0], row[1], row[2]
            # Skip multi-segment rows
            if "|" in left or "|" in right:
                continue
            # Strip trailing Ns from right (TSD trim padding)
            right_clean = right.rstrip("Nn")
            left_clean = left.rstrip("Nn").lstrip("Nn")
            # Reject if any internal Ns remain
            if "N" in left_clean.upper() or "N" in right_clean.upper():
                continue
            concat = (left_clean + right_clean).upper()
            if len(concat) < min_len:
                continue
            if concat in seen:
                continue
            seen.add(concat)
            seqs.append((acc, concat))
    return seqs

# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def run_palindrome_analysis(
    input_tsv, mge_name, output_dir,
    min_arm=3, max_spacer=2, max_mismatches=0,
    window_start=None, window_end=None,
    n_shuffles=100, seed=None, min_seq_len=20,
    centered=False,
):
    """Run palindrome detection + dinucleotide-shuffled null test.

    If centered=True, only the midpoint of the concatenated left+right
    sequence (i.e. the left/right flank boundary) is tested as the
    palindrome centre, for both observed sequences and all shuffled nulls.
    """
    logger = logging.getLogger("palindrome")

    if seed is not None:
        random.seed(seed)

    entries = load_sequences(input_tsv, min_len=min_seq_len)
    if not entries:
        logger.error(f"No valid sequences loaded from {input_tsv}")
        return

    seqs = [s for _, s in entries]
    total_len = len(seqs[0])
    # Sanity check: all same length?
    if not all(len(s) == total_len for s in seqs):
        logger.warning("Sequences have variable lengths; using max length for plots.")
        total_len = max(len(s) for s in seqs)

    # Window defaults to full sequence
    ws = 0 if window_start is None else max(0, window_start)
    we = total_len if window_end is None else min(total_len, window_end)

    # Fixed centre overrides window-based iteration when --centered is set.
    # The centre is the left/right flank boundary -- the midpoint of the
    # concatenated sequence -- not a hardcoded position, so this works
    # regardless of flank length.
    fixed_centre = total_len // 2 if centered else None

    logger.info(
        f"Loaded {len(seqs)} unique sequences (length {total_len} bp)"
    )
    logger.info(
        f"Pattern: arm >= {min_arm}, spacer <= {max_spacer}, "
        f"mismatches <= {max_mismatches}, window [{ws}, {we})"
        + (f", fixed centre = {fixed_centre}" if centered else "")
    )

   
    # --- Observed ---
    observed_hits = 0
    obs_positions = []
    obs_arm_lens = []
    obs_spacers = []
    for acc, s in entries:
        res = find_longest_inverted_repeat(
            s, max_spacer=max_spacer, min_arm=min_arm,
            max_mismatches=max_mismatches,
            search_start=ws, search_end=we,
            fixed_centre=fixed_centre,
        )
        if res["arm_len"] >= min_arm:
            observed_hits += 1
            obs_positions.append((
                res["left_start"], res["right_end"],
                res["arm_len"], res["spacer_size"],
            ))
            obs_arm_lens.append(res["arm_len"])
            obs_spacers.append(res["spacer_size"])

    logger.info(
        f"Observed: {observed_hits}/{len(seqs)} "
        f"({100 * observed_hits / len(seqs):.1f}%) hit the pattern"
    )

    # --- Null ---
    logger.info(f"Running {n_shuffles} dinucleotide shuffles...")
    null_hits = []
    null_positions = []
    null_arm_lens = []       # arm lens of hits only (for geometry stats)
    null_spacers = []
    null_arm_totals = []     # per-shuffle sum of arm lens across ALL sequences (non-hits = 0)
    pos_sample_cap = min(200, len(seqs))

    for i in range(n_shuffles):
        count = 0
        shuffle_arm_total = 0
        for idx, s in enumerate(seqs):
            try:
                shuffled = dinucleotide_shuffle(s)
            except RuntimeError:
                continue
            res = find_longest_inverted_repeat(
                shuffled, max_spacer=max_spacer, min_arm=min_arm,
                max_mismatches=max_mismatches,
                search_start=ws, search_end=we,
                fixed_centre=fixed_centre,
            )
            if res["arm_len"] >= min_arm:
                count += 1
                shuffle_arm_total += res["arm_len"]
                null_arm_lens.append(res["arm_len"])
                null_spacers.append(res["spacer_size"])
                if idx < pos_sample_cap:
                    null_positions.append((
                        res["left_start"], res["right_end"],
                        res["arm_len"], res["spacer_size"],
                    ))
        null_hits.append(count)
        null_arm_totals.append(shuffle_arm_total)
        if (i + 1) % 10 == 0:
            logger.debug(f"  {i+1}/{n_shuffles} shuffles done")

    # --- Stats ---
    null_mean = sum(null_hits) / len(null_hits) if null_hits else 0.0
    null_min = min(null_hits) if null_hits else 0
    null_max = max(null_hits) if null_hits else 0
    if len(null_hits) > 1:
        null_var = sum((x - null_mean) ** 2 for x in null_hits) / (len(null_hits) - 1)
        null_sd = null_var ** 0.5
    else:
        null_sd = 0.0
    p_val = (
        sum(n >= observed_hits for n in null_hits) / n_shuffles
        if n_shuffles > 0 else 1.0
    )
    # Hit-count enrichment (how many more sequences hit the threshold)
    hit_enrichment = observed_hits / null_mean if null_mean > 0 else float("inf")

    # Arm-length-weighted enrichment: mean arm len per sequence (non-hits = 0).
    # Reflects both frequency and strength, so a longer null arm deflates the score.
    obs_arm_total = sum(obs_arm_lens)
    obs_arm_mean_all = obs_arm_total / len(seqs)
    null_arm_mean_all = (sum(null_arm_totals) / len(null_arm_totals)) / len(seqs) if null_arm_totals else 0.0
    arm_enrichment = obs_arm_mean_all / null_arm_mean_all if null_arm_mean_all > 0 else float("inf")

    # Geometry stats
    def _mean(xs):
        return sum(xs) / len(xs) if xs else 0.0

    def _median(xs):
        if not xs:
            return 0.0
        s = sorted(xs)
        m = len(s) // 2
        return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2

    def _mode(xs):
        if not xs:
            return None
        counts = {}
        for x in xs:
            counts[x] = counts.get(x, 0) + 1
        return max(counts, key=counts.get)

    obs_arm_mean = _mean(obs_arm_lens)
    obs_arm_med = _median(obs_arm_lens)
    obs_arm_mode = _mode(obs_arm_lens)
    obs_sp_mean = _mean(obs_spacers)
    obs_sp_med = _median(obs_spacers)
    obs_sp_mode = _mode(obs_spacers)

    null_arm_mean = _mean(null_arm_lens)
    null_arm_med = _median(null_arm_lens)
    null_arm_mode = _mode(null_arm_lens)
    null_sp_mean = _mean(null_spacers)
    null_sp_med = _median(null_spacers)
    null_sp_mode = _mode(null_spacers)

    # Arm length breakdown (observed)
    obs_arm_dist = {}
    for a in obs_arm_lens:
        obs_arm_dist[a] = obs_arm_dist.get(a, 0) + 1
    null_arm_dist = {}
    for a in null_arm_lens:
        null_arm_dist[a] = null_arm_dist.get(a, 0) + 1

    summary_lines = [
        f"Palindrome / Inverted-Repeat Analysis: {mge_name}",
        "=" * 55,
        f"Input:                {input_tsv}",
        f"Sequences tested:     {len(seqs)} (unique concatenated L+R)",
        f"Sequence length:      {total_len} bp",
        f"Search window:        [{ws}, {we})",
        f"Centred mode:         {f'yes (centre fixed at index {fixed_centre})' if centered else 'no'}",
        f"Pattern minimum:      arm >= {min_arm}, spacer <= {max_spacer}, "
        f"mismatches <= {max_mismatches}",
        f"Null:                 {n_shuffles} dinucleotide shuffles",
        "",
        f"Observed hits:        {observed_hits}/{len(seqs)} "
        f"({100 * observed_hits / len(seqs):.1f}%)",
        f"Null mean hits:       {null_mean:.2f}/{len(seqs)} "
        f"({100 * null_mean / len(seqs):.1f}%)",
        f"Null SD:              {null_sd:.2f}",
        f"Null range:           {null_min}-{null_max}",
        f"Enrichment (hit rate):{hit_enrichment:>10.2f}x"
        if null_mean > 0 else "Enrichment (hit rate): inf",
        f"Enrichment (arm len): {arm_enrichment:>9.2f}x"
        if null_arm_mean_all > 0 else "Enrichment (arm len): inf",
        f"p-value (one-sided):  {p_val:.4f}"
        + (f"  (i.e. p < {1 / n_shuffles})" if p_val == 0 else ""),
        "",
        "Hit geometry (best IR per hit sequence):",
        f"{'':<22}{'observed':>12}{'null':>12}",
        f"{'  arm length mean':<22}{obs_arm_mean:>12.2f}{null_arm_mean:>12.2f}",
        f"{'  arm length median':<22}{obs_arm_med:>12.1f}{null_arm_med:>12.1f}",
        f"{'  arm length mode':<22}{str(obs_arm_mode):>12}{str(null_arm_mode):>12}",
        f"{'  spacer mean':<22}{obs_sp_mean:>12.2f}{null_sp_mean:>12.2f}",
        f"{'  spacer median':<22}{obs_sp_med:>12.1f}{null_sp_med:>12.1f}",
        f"{'  spacer mode':<22}{str(obs_sp_mode):>12}{str(null_sp_mode):>12}",
        "",
        "Arm length distribution (% of hits):",
        f"{'arm_len':>8}{'observed':>12}{'null':>12}",
    ]

    all_arms = sorted(set(list(obs_arm_dist.keys()) + list(null_arm_dist.keys())))
    for a in all_arms:
        obs_pct = 100 * obs_arm_dist.get(a, 0) / len(obs_arm_lens) if obs_arm_lens else 0
        null_pct = 100 * null_arm_dist.get(a, 0) / len(null_arm_lens) if null_arm_lens else 0
        summary_lines.append(f"{a:>8}{obs_pct:>11.1f}%{null_pct:>11.1f}%")

    logger.info(
        f"Palindrome result: {observed_hits}/{len(seqs)} ({100 * observed_hits / len(seqs):.1f}%) observed "
        f"vs null mean {null_mean:.2f} ({100 * null_mean / len(seqs):.1f}%); "
        f"enrichment {hit_enrichment:.2f}x (hit rate), {arm_enrichment:.2f}x (arm length); "
        f"p = {p_val:.4f} {_sig_stars(p_val)} (full breakdown written to file)"
    )
    for line in summary_lines:
        logger.debug(line)

    out_dir = Path(output_dir)
    files_dir = out_dir / "files"
    graphs_dir = out_dir / "graphs"
    files_dir.mkdir(parents=True, exist_ok=True)
    graphs_dir.mkdir(parents=True, exist_ok=True)

    summary_path = files_dir / f"{mge_name}_PALINDROME_SUMMARY.txt"
    try:
        with open(summary_path, "w") as f:
            f.write("\n".join(summary_lines) + "\n")
        logger.info(f"Summary written to {summary_path}")
    except Exception as e:
        logger.exception(f"Failed to write summary: {e}")

    # --- Plots ---
    try:
        null_dist_path = graphs_dir / f"{mge_name}_PALINDROME_NULL_DIST.png"
        plot_null_distribution(null_hits, observed_hits, p_val, mge_name, null_dist_path,
                                enrichment=hit_enrichment)
        logger.info(f"Null distribution plot written to {null_dist_path}")
    except Exception as e:
        logger.exception(f"Failed to write null distribution plot: {e}")

    try:
        arm_dist_path = graphs_dir / f"{mge_name}_PALINDROME_ARM_LEN_DIST.png"
        plot_arm_length_distribution(obs_arm_lens, null_arm_lens, n_shuffles, mge_name, arm_dist_path,
                                      p_val=p_val, enrichment=arm_enrichment)
        logger.info(f"Arm length distribution plot written to {arm_dist_path}")
    except Exception as e:
        logger.exception(f"Failed to write arm length distribution plot: {e}")

    logger.info("Palindrome analysis complete.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
    p = argparse.ArgumentParser(
        description="Detect inverted repeats in flank TSV and compare to dinucleotide-shuffled null."
    )
    p.add_argument("--input", "-i", required=True,
                   help="Trimmed flanks TSV: acc_name<TAB>left<TAB>right")
    p.add_argument("--mge_name", "-m", required=True,
                   help="MGE name used in output filenames")
    p.add_argument("--outdir", "-o", default=".",
                   help="Output directory (default: current)")
    p.add_argument("--min_arm", type=int, default=3,
                   help="Minimum arm length (default: 3)")
    p.add_argument("--max_spacer", type=int, default=2,
                   help="Maximum spacer size (default: 2)")
    p.add_argument("--max_mismatches", type=int, default=0,
                   help="Maximum mismatches per arm (default: 0)")
    p.add_argument("--window_start", type=int, default=None,
                   help="Window start (inclusive); default: 0")
    p.add_argument("--window_end", type=int, default=None,
                   help="Window end (exclusive); default: full sequence")
    p.add_argument("--n_shuffles", type=int, default=100,
                   help="Number of null shuffles (default: 100)")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed (optional)")
    p.add_argument("--min_seq_len", type=int, default=20,
                   help="Skip concatenated sequences shorter than this (default: 20)")
    p.add_argument("--centered", action="store_true",
                   help="Only test the left/right flank boundary as the palindrome centre (observed and null)")
    return p.parse_args()


def main():
    args = _parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        stream=sys.stdout,
    )
    run_palindrome_analysis(
        input_tsv=args.input,
        mge_name=args.mge_name,
        output_dir=args.outdir,
        min_arm=args.min_arm,
        max_spacer=args.max_spacer,
        max_mismatches=args.max_mismatches,
        window_start=args.window_start,
        window_end=args.window_end,
        n_shuffles=args.n_shuffles,
        seed=args.seed,
        min_seq_len=args.min_seq_len,
        centered=args.centered,
    )


if __name__ == "__main__":
    main()