"""
TSD (Target Site Duplication) Analysis Module.

Detects TSDs from paired flanking sequences, splits accessions into TSD /
no-TSD / trimmed groups, writes summary files, tests TSD enrichment against
a shuffled-pairing null, and generates plots for both the TSD hit rate and
the nucleotide distribution.

A TSD is a sequence duplicated at the end of the left flank and the start
of the right flank:
    Left:  XXXXXXAGGGTCA
    Right: AGGGTCAXXXXXX
"""

import logging
import random
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rc(seq):
    return seq.translate(_RC_TABLE)[::-1]


def _strip_n(seq):
    return seq.replace("N", "").replace("n", "")


def _write_lines(path, lines, label):
    """Write lines to a file, logging success or failure."""
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        logger.info(f"{label} written to {path}")
    except Exception as e:
        logger.exception(f"Failed to write {label}: {e}")


# ---------------------------------------------------------------------------
# TSD detection
# ---------------------------------------------------------------------------

def _find_longest_tsd(left_seq, right_seq, min_tsd=3):
    """Longest exact match between the end of left_seq and the start of right_seq (Ns ignored)."""
    best = 0
    for length in range(min_tsd, min(len(left_seq), len(right_seq)) + 1):
        l, r = _strip_n(left_seq[-length:]), _strip_n(right_seq[:length])
        if l and l == r:
            best = length
    return best


def _best_tsd_for_pairs(left_raw, right_raw, min_tsd=3):
    """Longest TSD across all left|right segment pairs for one accession. Returns (len, seq)."""
    left_segs, right_segs = left_raw.split("|"), right_raw.split("|")
    best_len, best_seq = 0, ""
    for l_seg, r_seg in zip(left_segs, right_segs):
        if l_seg == "none" or r_seg == "none":
            continue
        tsd_len = _find_longest_tsd(l_seg, r_seg, min_tsd=min_tsd)
        if tsd_len > best_len:
            best_len, best_seq = tsd_len, l_seg[-tsd_len:]
    return best_len, best_seq, left_segs, right_segs


def find_tsds(left_seqs, right_seqs, acc_names, min_tsd=3):
    """Find the longest TSD per accession, pairing left/right segments by index.

    Multi-flank entries ("a|b|c") are split on "|" and paired positionally.
    Returns a list of dicts: acc_name, pairs, tsd_len, tsd_seq.
    """
    results = []
    for acc, left_raw, right_raw in zip(acc_names, left_seqs, right_seqs):
        best_len, best_seq, left_segs, right_segs = _best_tsd_for_pairs(left_raw, right_raw, min_tsd)
        if len(left_segs) != len(right_segs):
            logger.warning(
                f"{acc}: {len(left_segs)} left vs {len(right_segs)} right segments — "
                f"pairing up to the shorter count"
            )
        results.append({
            "acc_name": acc,
            "pairs": list(zip(left_segs, right_segs)),
            "tsd_len": best_len,
            "tsd_seq": best_seq,
        })

    hits = sum(1 for r in results if r["tsd_len"] > 0)
    logger.info(f"TSD detection: {hits}/{len(results)} accessions have a TSD")
    return results


# ---------------------------------------------------------------------------
# Null model: break true left/right pairing by shuffling, recheck for TSDs
# ---------------------------------------------------------------------------

def shuffled_null_tsd_lengths(left_seqs, right_seqs, min_tsd=3, n_shuffles=100, seed=None):
    """Shuffle right_seqs against left_seqs to break real pairing, repeat n_shuffles times.

    Returns (null_hit_counts, null_tsd_lengths):
      null_hit_counts: hits per shuffle (for the significance test)
      null_tsd_lengths: TSD length of every hit across all shuffles (for the length plot)
    """
    if seed is not None:
        random.seed(seed)

    null_hit_counts, null_tsd_lengths = [], []
    right_pool = list(right_seqs)

    for _ in range(n_shuffles):
        shuffled_right = right_pool.copy()
        random.shuffle(shuffled_right)
        hits = 0
        for left_raw, right_raw in zip(left_seqs, shuffled_right):
            tsd_len, _, _, _ = _best_tsd_for_pairs(left_raw, right_raw, min_tsd)
            if tsd_len > 0:
                hits += 1
                null_tsd_lengths.append(tsd_len)
        null_hit_counts.append(hits)

    return null_hit_counts, null_tsd_lengths


def tsd_significance(tsd_results, left_seqs, right_seqs, min_tsd=3, n_shuffles=100, seed=None):
    """Compare observed TSD hit rate to a shuffled-pairing null. Returns a stats dict."""
    observed_hits = sum(1 for r in tsd_results if r["tsd_len"] > 0)
    obs_lengths = [r["tsd_len"] for r in tsd_results if r["tsd_len"] > 0]

    null_hit_counts, null_lengths = shuffled_null_tsd_lengths(
        left_seqs, right_seqs, min_tsd=min_tsd, n_shuffles=n_shuffles, seed=seed
    )
    null_mean = sum(null_hit_counts) / len(null_hit_counts) if null_hit_counts else 0.0
    p_val = sum(n >= observed_hits for n in null_hit_counts) / n_shuffles if n_shuffles else 1.0
    enrichment = observed_hits / null_mean if null_mean > 0 else float("inf")

    return {
        "observed_hits": observed_hits,
        "obs_lengths": obs_lengths,
        "null_hit_counts": null_hit_counts,
        "null_lengths": null_lengths,
        "null_mean": null_mean,
        "enrichment": enrichment,
        "p_value": p_val,
        "n_shuffles": n_shuffles,
    }


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


def plot_tsd_null_distribution(stats, mge_name, output_path):
    """Histogram of per-shuffle TSD hit counts, with the observed value marked."""
    null_hits, observed = stats["null_hit_counts"], stats["observed_hits"]
    p_val, enrichment = stats["p_value"], stats["enrichment"]
    stars = _sig_stars(p_val)

    fig, ax = plt.subplots(figsize=(8, 5))
    lo, hi = min(null_hits + [observed]), max(null_hits + [observed])
    ax.hist(null_hits, bins=range(lo, hi + 2), color="tab:orange",
            edgecolor="white", align="left", label="Null (shuffled pairing)")
    ax.axvline(observed, color="tab:blue", linewidth=2.5, label=f"Observed ({observed})")
    ax.set_xlabel("Accessions with a TSD")
    ax.set_ylabel("Number of shuffles")
    ax.set_title(f"{mge_name}: observed vs. null TSD hits (p = {p_val:.4f} {stars})")
    ax.text(0.98, 0.95, f"Enrichment: {enrichment:.2f}x\np = {p_val:.4f} {stars}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_tsd_length_distribution(stats, mge_name, output_path):
    """Bar chart of TSD length distribution, observed vs. null (averaged per shuffle).

    The annotated p-value/enrichment are the overall hit-rate significance
    (from plot_tsd_null_distribution's test), not a per-length test -- this
    chart only breaks that overall result down by TSD length.
    """
    obs_lengths, null_lengths, n_shuffles = stats["obs_lengths"], stats["null_lengths"], stats["n_shuffles"]
    all_lengths = sorted(set(obs_lengths) | set(null_lengths))
    if not all_lengths:
        return
    obs_counts = [obs_lengths.count(l) for l in all_lengths]
    null_counts = [null_lengths.count(l) / n_shuffles for l in all_lengths]

    p_val, enrichment = stats["p_value"], stats["enrichment"]
    stars = _sig_stars(p_val)

    x = range(len(all_lengths))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([i - width / 2 for i in x], obs_counts, width, label="Observed", color="tab:blue")
    ax.bar([i + width / 2 for i in x], null_counts, width,
           label="Null (avg per shuffle)", color="tab:orange")
    ax.set_xticks(list(x))
    ax.set_xticklabels(all_lengths)
    ax.set_xlabel("TSD length (bp)")
    ax.set_ylabel("Number of hits")
    ax.set_title(f"{mge_name}: TSD length, observed vs. null")
    ax.text(0.98, 0.95, f"Overall enrichment: {enrichment:.2f}x\np = {p_val:.4f} {stars}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.legend(loc="upper left")
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


# ---------------------------------------------------------------------------
# Splitting and file writing
# ---------------------------------------------------------------------------

def build_split_entries(tsd_results, flank_len):
    """Split TSD results into three deduplicated groups.

    Returns (tsd_entries, no_tsd_entries, trimmed_entries), each a list of
    dicts with keys: acc_name, left, right.

    Trimmed entries have the TSD removed from the right flank, N-padded on
    the right to preserve length.
    """
    tsd_entries, no_tsd_entries, trimmed_entries = [], [], []
    seen_tsd, seen_no_tsd, seen_trimmed = set(), set(), set()

    for result in tsd_results:
        acc, tsd_len = result["acc_name"], result["tsd_len"]

        left_rep, right_rep = next(
            ((l, r) for l, r in result["pairs"] if l != "none" and r != "none"),
            (None, None)
        )
        if left_rep is None:
            logger.debug(f"No valid pairs for {acc}, skipping")
            continue

        key = (left_rep, right_rep)
        if tsd_len == 0:
            if key not in seen_no_tsd:
                seen_no_tsd.add(key)
                no_tsd_entries.append({"acc_name": acc, "left": left_rep, "right": right_rep})
            continue

        if key not in seen_tsd:
            seen_tsd.add(key)
            tsd_entries.append({"acc_name": acc, "left": left_rep, "right": right_rep})

        padded_right = right_rep[tsd_len:] + "N" * tsd_len
        if len(padded_right) != flank_len:
            logger.warning(
                f"Padded right flank length error for {acc}: "
                f"expected {flank_len}, got {len(padded_right)}"
            )
        trimmed_key = (left_rep, padded_right)
        if trimmed_key not in seen_trimmed:
            seen_trimmed.add(trimmed_key)
            trimmed_entries.append({"acc_name": acc, "left": left_rep, "right": padded_right})

    logger.info(
        f"Split: {len(tsd_entries)} TSD, {len(no_tsd_entries)} no-TSD, "
        f"{len(trimmed_entries)} trimmed (all unique)"
    )
    return tsd_entries, no_tsd_entries, trimmed_entries


def _write_summary_file(path, entries, flank_len):
    """Write TETyper-format summary file (acc\\tleft\\tright)."""
    for e in entries:
        for side, seq in [("left", e["left"]), ("right", e["right"])]:
            if len(seq) != flank_len:
                logger.warning(
                    f"{side} flank length mismatch for {e['acc_name']}: "
                    f"expected {flank_len}, got {len(seq)}"
                )
    lines = [f"{e['acc_name']}\t{e['left']}\t{e['right']}" for e in entries]
    _write_lines(path, lines, f"{Path(path).name}")


# ---------------------------------------------------------------------------
# TSD summary
# ---------------------------------------------------------------------------

def write_tsd_summary(tsd_results, mge_name, output_dir):
    """Log and save a TSD detection summary (counts, length breakdown, per-accession detail)."""
    total = len(tsd_results)
    tsd_count = sum(1 for r in tsd_results if r["tsd_len"] > 0)
    length_counts = Counter(r["tsd_len"] for r in tsd_results if r["tsd_len"] > 0)

    header = [
        f"TSD Summary for {mge_name}",
        "=" * 40,
        f"Total accessions: {total}",
        f"With TSD:         {tsd_count} ({tsd_count/total:.1%})",
        f"Without TSD:      {total - tsd_count} ({(total - tsd_count)/total:.1%})",
        "",
        "TSD length breakdown:",
    ] + [
        f"  Length {l:>3}: {c:>6} accessions ({c/total:.1%})"
        for l, c in sorted(length_counts.items())
    ]
    detail = [
        "",
        "Per-accession details:",
        f"{'acc_name':<30}\t{'tsd_len':<10}\ttsd_seq",
    ] + [
        f"{r['acc_name']:<30}\t{r['tsd_len']:<10}\t{r['tsd_seq']}"
        for r in tsd_results
    ]

    logger.info(
        f"TSD summary: {tsd_count}/{total} ({tsd_count/total:.1%}) with TSD, "
        f"{total - tsd_count} ({(total - tsd_count)/total:.1%}) without "
        "(full length breakdown and per-accession detail in file)"
    )
    for line in header:
        logger.debug(line)

    _write_lines(Path(output_dir) / "files" / f"{mge_name}_TSD_SUMMARY.txt", header + detail, "TSD summary")


def run_tsd_analysis(
    left_seqs, right_seqs, acc_names, flank_len, mge_name, output_dir,
    min_tsd=3, n_shuffles=100, seed=None, test_significance=False,
):
    """TSD analysis pipeline: detect, split, write files, and (optionally,
    since it's the expensive part -- O(n_seqs) per shuffle) test the hit
    rate against a shuffled-pairing null and plot it.
    """
    files_dir = Path(output_dir) / "files"
    graphs_dir = Path(output_dir) / "graphs"

    tsd_results = find_tsds(left_seqs, right_seqs, acc_names, min_tsd=min_tsd)
    write_tsd_summary(tsd_results, mge_name, output_dir)

    stats = {
        'observed_hits': sum(1 for r in tsd_results if r['tsd_len'] > 0),
        'total': len(tsd_results),
        'length_counts': Counter(r['tsd_len'] for r in tsd_results if r['tsd_len'] > 0),
    }

    if test_significance:
        sig_stats = tsd_significance(
            tsd_results, left_seqs, right_seqs,
            min_tsd=min_tsd, n_shuffles=n_shuffles, seed=seed,
        )
        stats.update(sig_stats)
        logger.info(
            f"TSD significance: observed {stats['observed_hits']}, "
            f"null mean {stats['null_mean']:.2f}, enrichment {stats['enrichment']:.2f}x, "
            f"p = {stats['p_value']:.4f} {_sig_stars(stats['p_value'])}"
            + (f" (i.e. p < {1 / n_shuffles})" if stats["p_value"] == 0 else "")
        )
        try:
            graphs_dir.mkdir(parents=True, exist_ok=True)
            plot_tsd_null_distribution(stats, mge_name, graphs_dir / f"{mge_name}_TSD_NULL_DIST.png")
            plot_tsd_length_distribution(stats, mge_name, graphs_dir / f"{mge_name}_TSD_LENGTH_DIST.png")
            logger.info(f"TSD significance plots written to {graphs_dir}")
        except Exception as e:
            logger.exception(f"Failed to write TSD significance plots: {e}")
    else:
        logger.info(
            f"TSD detection: {stats['observed_hits']}/{stats['total']} "
            f"({stats['observed_hits'] / stats['total']:.1%}) accessions have a TSD "
            "(significance test skipped)"
        )

    tsd_entries, no_tsd_entries, trimmed_entries = build_split_entries(tsd_results, flank_len)

    _write_summary_file(files_dir / f"{mge_name}_TSD.txt",         tsd_entries,     flank_len)
    _write_summary_file(files_dir / f"{mge_name}_no_TSD.txt",      no_tsd_entries,  flank_len)

    logger.info("TSD analysis complete.")
    return stats