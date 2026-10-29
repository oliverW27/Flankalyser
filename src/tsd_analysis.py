"""
TSD (Target Site Duplication) Analysis Module.

Detects TSDs in each unique left/right flank pair and writes a summary file.
Optionally tests the TSD hit rate against the same flanks with the boundary
shifted away from the insertion site, and plots it by TSD length.

A TSD is a sequence duplicated at the end of the left flank and the start
of the right flank:
    Left:  XXXXXXAGGGTCA
    Right: AGGGTCAXXXXXX
"""

import logging
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from scipy.stats import binom

logger = logging.getLogger(__name__)

MAX_TSD = 25

# Background rate comes from re-testing each pair with the boundary moved this many bases
# into the left (+) or right (-) flank. A real TSD of length T still gives a shorter exact
# match when shifted by up to T - min_tsd bases, so shifts start well beyond typical TSDs.
NULL_SHIFTS = [*range(-25, -9), *range(10, 26)]


def _write_lines(path, lines, label):
    """Write lines to a file, logging success or failure."""
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        logger.info(f"{label} written to {path}")
    except Exception as e:
        logger.exception(f"Failed to write {label}: {e}")


def unique_pairs(left_data, right_data, acc_names, flank_len):
    """Unique (left, right) flank pairs.

    Each insertion site counts once however many genomes it appears in, as in
    the rest of the pipeline. Pipe-separated segments are paired by position;
    'none' and wrong-length segments are skipped.

    Returns:
        list of (left, right)
    """
    pairs = {}  # dict keeps first-seen order
    for acc, left_field, right_field in zip(acc_names, left_data, right_data):
        left_segs, right_segs = left_field.upper().split("|"), right_field.upper().split("|")
        if len(left_segs) != len(right_segs):
            logger.warning(
                f"{acc}: {len(left_segs)} left vs {len(right_segs)} right segments -- "
                f"pairing up to the shorter count"
            )
        for left, right in zip(left_segs, right_segs):
            if len(left) == flank_len and len(right) == flank_len:
                pairs[(left, right)] = None
    return list(pairs)


def _find_longest_tsd(left_seq, right_seq, min_tsd=3):
    """Longest exact duplication between the end of left_seq and the start of right_seq.

    A TSD may not contain N (an N never counts as matching anything);
    Ns elsewhere in the flanks don't matter.
    """
    best = 0
    for length in range(min_tsd, min(MAX_TSD, len(left_seq), len(right_seq)) + 1):
        tsd = left_seq[-length:]
        if tsd == right_seq[:length] and "N" not in tsd:
            best = length
    return best


def find_tsds(pairs, min_tsd=3):
    """Find the longest TSD in each unique pair.

    Returns a list of dicts: left, right, tsd_len.
    """
    return [{"left": left, "right": right, "tsd_len": _find_longest_tsd(left, right, min_tsd)}
            for left, right in pairs]


def shifted_null_tsd_lengths(tsd_results, min_tsd=3):
    """Look for a TSD at each of NULL_SHIFTS shifted boundaries in every pair.

    Returns (null_hit_counts, null_tsd_lengths):
      null_hit_counts: hits per shift (for the significance test)
      null_tsd_lengths: TSD length of every hit across all shifts (for the per-length background)
    """
    null_hit_counts, null_tsd_lengths = [], []
    for shift in NULL_SHIFTS:
        hits = []
        for r in tsd_results:
            left, right = r["left"], r["right"]
            if shift > 0:
                left = left[:-shift]
            else:
                right = right[-shift:]
            tsd_len = _find_longest_tsd(left, right, min_tsd)
            if tsd_len:
                hits.append(tsd_len)
        null_hit_counts.append(len(hits))
        null_tsd_lengths.extend(hits)
    return null_hit_counts, null_tsd_lengths


def tsd_significance(tsd_results, min_tsd=3):
    """Compare the TSD hit rate at the real boundary to the rate at shifted boundaries.

    p-value is a one-sided binomial test of the observed hits against the
    background rate pooled over all shifts.
    """
    observed_hits = sum(1 for r in tsd_results if r["tsd_len"] > 0)
    n = len(tsd_results)

    null_hit_counts, null_lengths = shifted_null_tsd_lengths(tsd_results, min_tsd=min_tsd)
    null_mean = sum(null_hit_counts) / len(null_hit_counts)
    return {
        "null_mean": null_mean,
        "enrichment": observed_hits / null_mean if null_mean > 0 else float("inf"),
        "p_value": float(binom.sf(observed_hits - 1, n, null_mean / n)) if null_mean > 0 else 0.0,
        # average number of sites per shifted boundary with a TSD of each length
        "null_length_counts": {l: c / len(NULL_SHIFTS) for l, c in Counter(null_lengths).items()},
    }


def _fmt_p(p):
    """p-value for display; a binomial p this small underflows to exactly 0."""
    return "< 1e-300" if p == 0 else f"= {p:.2g}"


def _sig_stars(p):
    """Conventional significance stars for a one-sided p-value."""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_tsd_length_distribution(stats, mge_name, output_path):
    """Bar chart of TSD length distribution, observed vs. null (averaged per shifted boundary).

    The annotated p-value/enrichment are the overall hit-rate significance,
    not a per-length test -- this chart only breaks that overall result down
    by TSD length.
    """
    obs_by_len, null_by_len = stats["length_counts"], stats["null_length_counts"]
    all_lengths = sorted(set(obs_by_len) | set(null_by_len))
    if not all_lengths:
        return
    obs_counts = [obs_by_len.get(l, 0) for l in all_lengths]
    null_counts = [null_by_len.get(l, 0) for l in all_lengths]

    p_val, enrichment = stats["p_value"], stats["enrichment"]
    stars = _sig_stars(p_val)

    x = range(len(all_lengths))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([i - width / 2 for i in x], obs_counts, width, label="Observed", color="tab:blue")
    ax.bar([i + width / 2 for i in x], null_counts, width,
           label="Null (avg per shifted boundary)", color="tab:orange")
    ax.set_xticks(list(x))
    ax.set_xticklabels(all_lengths)
    ax.set_xlabel("TSD length (bp)")
    ax.set_ylabel("Number of hits")
    ax.set_title(f"{mge_name}: TSD length, observed vs. null")
    ax.text(0.98, 0.95, f"Overall enrichment: {enrichment:.2f}x\np {_fmt_p(p_val)} {stars}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def write_tsd_summary(stats, mge_name, output_dir):
    """Save the TSD summary: counts, significance (if tested), and length breakdown."""
    total, tsd_count = stats["total"], stats["observed_hits"]
    background = stats.get("null_length_counts")

    if background is None:
        significance = "Significance:           not tested (use --tsd_test_significance)"
    else:
        significance = (
            f"Significance:           background {stats['null_mean']:.1f} sites "
            f"(boundary shifted 10-25 bp), enrichment {stats['enrichment']:.2f}x, "
            f"p {_fmt_p(stats['p_value'])} {_sig_stars(stats['p_value'])}"
        )

    lengths = sorted(set(stats["length_counts"]) | set(background or {}))
    breakdown = []
    for l in lengths:
        c = stats["length_counts"].get(l, 0)
        line = f"  Length {l:>3}: {c:>6} sites ({c/total:.1%})"
        if background is not None:
            line += f"   background {background.get(l, 0):.1f}"
        breakdown.append(line)

    lines = [
        f"TSD Summary for {mge_name}",
        "=" * 40,
        f"Unique insertion sites: {total}",
        f"With TSD:               {tsd_count} ({tsd_count/total:.1%})",
        f"Without TSD:            {total - tsd_count} ({(total - tsd_count)/total:.1%})",
        significance,
        "",
        "TSD length breakdown" + (" (background = average sites per shifted boundary):"
                                  if background is not None else ":"),
    ] + breakdown

    _write_lines(Path(output_dir) / "files" / f"{mge_name}_TSD_SUMMARY.txt", lines, "TSD summary")


def run_tsd_analysis(
    left_data, right_data, acc_names, flank_len, mge_name, output_dir,
    min_tsd=3, test_significance=False,
):
    """TSD analysis pipeline: deduplicate, detect, optionally test the hit rate
    against shifted boundaries (and plot it by length), and write the summary.

    left_data / right_data / acc_names are the raw per-line fields from the summary file.
    """
    graphs_dir = Path(output_dir) / "graphs"

    pairs = unique_pairs(left_data, right_data, acc_names, flank_len)
    if not pairs:
        logger.warning("No usable left/right flank pairs -- skipping TSD analysis")
        return
    tsd_results = find_tsds(pairs, min_tsd=min_tsd)

    length_counts = Counter(r['tsd_len'] for r in tsd_results if r['tsd_len'] > 0)
    stats = {
        'observed_hits': sum(length_counts.values()),
        'total': len(tsd_results),
        'length_counts': length_counts,
    }
    logger.info(
        f"TSD detection: {stats['observed_hits']}/{stats['total']} "
        f"({stats['observed_hits'] / stats['total']:.1%}) unique insertion sites have a TSD"
    )

    if test_significance:
        stats.update(tsd_significance(tsd_results, min_tsd=min_tsd))
        logger.info(
            f"TSD significance: shifted-boundary background {stats['null_mean']:.1f}, "
            f"enrichment {stats['enrichment']:.2f}x, "
            f"p {_fmt_p(stats['p_value'])} {_sig_stars(stats['p_value'])}"
        )
        try:
            graphs_dir.mkdir(parents=True, exist_ok=True)
            plot_tsd_length_distribution(stats, mge_name, graphs_dir / f"{mge_name}_TSD_LENGTH_DIST.png")
            logger.info(f"TSD length plot written to {graphs_dir}")
        except Exception as e:
            logger.exception(f"Failed to write TSD length plot: {e}")

    write_tsd_summary(stats, mge_name, output_dir)
