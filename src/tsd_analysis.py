"""
TSD (Target Site Duplication) Analysis Module.

Detects TSDs from paired flanking sequences, splits accessions into TSD /
no-TSD / trimmed groups, writes summary files, generates nucleotide
distribution graphs, and checks trimmed sequences for REP-like stem-loops.

A TSD is a sequence duplicated at the end of the left flank and the start
of the right flank:
    Left:  XXXXXXAGGGTCA
    Right: AGGGTCAXXXXXX
"""

import logging
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.nucl_dist import count_nucl, nucl_dist_stat_sig, calc_nucl_props, graph, adjust_gc_content

logger = logging.getLogger(__name__)

_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rc(seq):
    return seq.translate(_RC_TABLE)[::-1]


def _strip_n(seq):
    return seq.replace("N", "").replace("n", "")


def _write_file(path, lines, label, logger):
    """Write lines to a file, logging success or failure."""
    try:
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        logger.info(f"{label} written to {path}")
    except Exception as e:
        logger.exception(f"Failed to write {label}: {e}")


# ---------------------------------------------------------------------------
# TSD detection
# ---------------------------------------------------------------------------

def _find_longest_tsd(left_seq, right_seq, min_tsd=3):
    """Return length of longest TSD between left suffix and right prefix (Ns ignored)."""
    best = 0
    for length in range(min_tsd, min(len(left_seq), len(right_seq)) + 1):
        l = _strip_n(left_seq[-length:])
        r = _strip_n(right_seq[:length])
        if l and r and l == r:
            best = length
    return best


def find_tsds(left_seqs, right_seqs, acc_names, flank_len, logger, min_tsd=3):
    """Find the longest TSD per accession, pairing L_i with R_i for multi-flank entries.

    Returns a list of dicts with keys: acc_name, pairs, tsd_len, tsd_seq.
    """
    results = []
    for acc, left_raw, right_raw in zip(acc_names, left_seqs, right_seqs):
        left_segs = left_raw.split("|")
        right_segs = right_raw.split("|")
        if len(left_segs) != len(right_segs):
            logger.warning(
                f"Unequal segment counts for {acc}: "
                f"{len(left_segs)} left vs {len(right_segs)} right — "
                f"pairing up to the shorter count"
            )
        best_len, best_seq = 0, ""
        for l_seg, r_seg in zip(left_segs, right_segs):
            if l_seg == "none" or r_seg == "none":
                continue
            tsd_len = _find_longest_tsd(l_seg, r_seg, min_tsd=min_tsd)
            if tsd_len > best_len:
                best_len = tsd_len
                best_seq = l_seg[-tsd_len:]
        results.append({
            "acc_name": acc,
            "pairs": list(zip(left_segs, right_segs)),
            "tsd_len": best_len,
            "tsd_seq": best_seq,
        })

    tsd_count = sum(1 for r in results if r["tsd_len"] > 0)
    logger.info(f"TSD detection complete: {tsd_count}/{len(results)} accessions have a TSD")
    return results


# ---------------------------------------------------------------------------
# Splitting and file writing
# ---------------------------------------------------------------------------

def build_split_entries(tsd_results, flank_len, logger):
    """Split TSD results into three deduplicated groups.

    Returns (tsd_entries, no_tsd_entries, trimmed_entries), each a list of
    dicts with keys: acc_name, left, right.

    Trimmed entries have the TSD removed from the right flank, N-padded on
    the right to preserve length.
    """
    tsd_entries, no_tsd_entries, trimmed_entries = [], [], []
    seen_tsd, seen_no_tsd, seen_trimmed = set(), set(), set()

    for result in tsd_results:
        acc = result["acc_name"]
        tsd_len = result["tsd_len"]

        # First valid pair is the representative
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
        else:
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


def _write_summary_file(path, entries, flank_len, logger):
    """Write TETyper-format summary file (acc\\tleft\\tright)."""
    try:
        with open(path, "w") as f:
            for e in entries:
                for side, seq in [("left", e["left"]), ("right", e["right"])]:
                    if len(seq) != flank_len:
                        logger.warning(
                            f"{side} flank length mismatch for {e['acc_name']}: "
                            f"expected {flank_len}, got {len(seq)}"
                        )
                f.write(f"{e['acc_name']}\t{e['left']}\t{e['right']}\n")
        logger.info(f"Wrote {len(entries)} entries to {path}")
    except Exception as e:
        logger.exception(f"Failed to write {path}: {e}")
        raise


# ---------------------------------------------------------------------------
# TSD summary
# ---------------------------------------------------------------------------

def write_tsd_summary(tsd_results, mge_name, output_dir, logger):
    """Log and save TSD detection summary (counts, length breakdown, per-accession detail)."""
    total = len(tsd_results)
    tsd_found = [r for r in tsd_results if r["tsd_len"] > 0]
    tsd_count = len(tsd_found)

    length_counts = {}
    for r in tsd_found:
        length_counts[r["tsd_len"]] = length_counts.get(r["tsd_len"], 0) + 1

    lines = [
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
    ] + [
        "",
        "Per-accession details:",
        f"{'acc_name':<30}\t{'tsd_len':<10}\ttsd_seq",
    ] + [
        f"{r['acc_name']:<30}\t{r['tsd_len']:<10}\t{r['tsd_seq']}"
        for r in tsd_results
    ]

    for line in lines[:8 + len(length_counts)]:
        logger.info(line)

    _write_file(
        Path(output_dir) / "files" / f"{mge_name}_TSD_summary.txt",
        lines, "TSD summary", logger
    )


# ---------------------------------------------------------------------------
# Nucleotide graph runner
# ---------------------------------------------------------------------------

def _run_nucl_graph(seqs_left, seqs_right, label, flank_len, output_dir,
                    adjust_proportions=False, adjust_gc=False, gc_content=None):
    """Run nucleotide distribution counting and graph generation for one group."""
    if not seqs_left and not seqs_right:
        logger.warning(f"No sequences for '{label}', skipping graph")
        return
    try:
        left_counts, right_counts = count_nucl(seqs_left, seqs_right, flank_len)
        left_stats = nucl_dist_stat_sig(left_counts, gc_content=gc_content)
        right_stats = nucl_dist_stat_sig(right_counts, gc_content=gc_content)

        if adjust_proportions:
            combined_props = calc_nucl_props(left_counts, right_counts, flank_len)
        else:
            default_prop = [0.25] * flank_len + ["nan"] + [0.25] * flank_len
            combined_props = {n: list(default_prop) for n in ["A", "C", "G", "T"]}

        if adjust_gc:
            combined_props = adjust_gc_content(combined_props, flank_len)

        graph(
            flank_len, label, combined_props, left_stats, right_stats,
            outdir=output_dir,
            total_left=len(seqs_left), unique_left=len(set(seqs_left)),
            total_right=len(seqs_right), unique_right=len(set(seqs_right)),
            gc_content=gc_content,
        )
        logger.info(f"Graph generated for '{label}'")
    except Exception as e:
        logger.exception(f"Graph failed for '{label}': {e}")


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def run_tsd_analysis(
    left_seqs, right_seqs, acc_names, flank_len, mge_name, output_dir, logger,
    adjust_proportions=False, adjust_gc=False, gc_content=None, min_tsd=3,
):
    """Full TSD analysis pipeline: detect, summarise, split, write files, graph."""
    files_dir = Path(output_dir) / "files"

    tsd_results = find_tsds(left_seqs, right_seqs, acc_names, flank_len, logger, min_tsd=min_tsd)
    write_tsd_summary(tsd_results, mge_name, output_dir, logger)

    tsd_entries, no_tsd_entries, trimmed_entries = build_split_entries(
        tsd_results, flank_len, logger
    )

    _write_summary_file(str(files_dir / f"{mge_name}_TSD.txt"),         tsd_entries,      flank_len, logger)
    _write_summary_file(str(files_dir / f"{mge_name}_no_TSD.txt"),      no_tsd_entries,   flank_len, logger)
    _write_summary_file(str(files_dir / f"{mge_name}_trimmed_TSD.txt"), trimmed_entries,  flank_len, logger)

    for label, entries in [
        (f"{mge_name}_TSD",         tsd_entries),
        (f"{mge_name}_NO-TSD",      no_tsd_entries),
        (f"{mge_name}_TRIMMED-TSD", trimmed_entries),
    ]:
        _run_nucl_graph(
            [e["left"] for e in entries], [e["right"] for e in entries],
            label, flank_len, output_dir,
            adjust_proportions=adjust_proportions,
            adjust_gc=adjust_gc,
            gc_content=gc_content,
        )

    logger.info("TSD analysis complete.")