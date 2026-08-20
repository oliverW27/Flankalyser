# Tests for src/tsd_analysis.py
import logging

import pytest

from src.tsd_analysis import (
    _best_tsd_for_pairs,
    _find_longest_tsd,
    build_split_entries,
    find_tsds,
    run_tsd_analysis,
    shuffled_null_tsd_lengths,
    tsd_significance,
    write_tsd_summary,
)


def _make_flank_pair(flank_len, tsd_seq, left_filler="A", right_filler="C"):
    """Build one (left, right) flank pair with tsd_seq planted at the boundary:
    left ends with it, right starts with it. Fillers are single, distinct
    repeated bases so no *other* shared run can appear anywhere else in the
    pair by coincidence.
    """
    k = len(tsd_seq)
    left = left_filler * (flank_len - k) + tsd_seq
    right = tsd_seq + right_filler * (flank_len - k)
    return left, right


def test_find_longest_tsd_exact_match():
    """A right flank starting with the same k bases that end the left flank
    should return the correct (longest) TSD length."""
    left, right = _make_flank_pair(flank_len=20, tsd_seq="GCATGCA")
    assert _find_longest_tsd(left, right, min_tsd=3) == 7


def test_find_longest_tsd_below_min_tsd_returns_zero():
    """A shared suffix/prefix shorter than min_tsd should return 0, not the
    short match."""
    left, right = _make_flank_pair(flank_len=20, tsd_seq="GC")
    assert _find_longest_tsd(left, right, min_tsd=5) == 0
    assert _find_longest_tsd(left, right, min_tsd=2) == 2


def test_find_longest_tsd_ignores_n_characters():
    """N's are stripped from each window before comparison. If the N sits at
    the same relative position in both windows, both windows shrink by one
    and the comparison still succeeds at the full length. (If the N only
    appears on one side, the two stripped windows end up different lengths
    and that comparison fails -- the match then falls back to a shorter,
    N-free window instead.)"""
    left, right = _make_flank_pair(flank_len=20, tsd_seq="GCANGCA")  # symmetric N
    assert _find_longest_tsd(left, right, min_tsd=3) == 7

    # N on only one side: the length-7 comparison fails (7 vs 6 chars after
    # stripping), but the match still recovers at the longest N-free window.
    left_asym = "A" * 13 + "GCATGCA"
    right_asym = "GCANGCA" + "C" * 13
    assert _find_longest_tsd(left_asym, right_asym, min_tsd=3) == 3


def test_best_tsd_for_pairs_skips_none_segments():
    """Segment pairs where either side is 'none' should be skipped entirely."""
    tsd_left, tsd_right = _make_flank_pair(flank_len=20, tsd_seq="GCATGCA")
    left_raw = f"none|{tsd_left}"
    right_raw = f"{tsd_right}|none"

    best_len, best_seq, left_segs, right_segs = _best_tsd_for_pairs(left_raw, right_raw, min_tsd=3)

    # Pair 0 is (none, tsd_right) and pair 1 is (tsd_left, none) -- both have
    # a 'none' side, so neither is ever compared.
    assert best_len == 0
    assert best_seq == ""
    assert left_segs == ["none", tsd_left]
    assert right_segs == [tsd_right, "none"]


def test_best_tsd_for_pairs_picks_longest_across_multiple_segments():
    """When multiple left|right segment pairs exist, the best (longest) TSD
    across all pairs should be selected, with the correct tsd_seq returned."""
    short_left, short_right = _make_flank_pair(flank_len=20, tsd_seq="GCA")     # 3bp TSD
    long_left, long_right = _make_flank_pair(flank_len=20, tsd_seq="GCATGCA")   # 7bp TSD

    left_raw = f"{short_left}|{long_left}"
    right_raw = f"{short_right}|{long_right}"

    best_len, best_seq, left_segs, right_segs = _best_tsd_for_pairs(left_raw, right_raw, min_tsd=3)

    assert best_len == 7
    assert best_seq == "GCATGCA"


def test_find_tsds_mismatched_segment_counts_logs_warning(caplog):
    """Unequal numbers of left vs right pipe-segments should log a warning
    and pair up to the shorter count (via zip) without crashing."""
    tsd_left, tsd_right = _make_flank_pair(flank_len=20, tsd_seq="GCATGCA")
    # 2 left segments, only 1 right segment
    left_raw = f"{tsd_left}|{tsd_left}"
    right_raw = tsd_right

    with caplog.at_level(logging.WARNING):
        results = find_tsds([left_raw], [right_raw], ["acc1"], min_tsd=3)

    assert any("2 left vs 1 right segments" in record.message for record in caplog.records)
    # Still pairs up to the shorter count and finds the TSD rather than crashing.
    assert results[0]["tsd_len"] == 7


def test_shuffled_null_tsd_lengths_reproducible_with_seed():
    """Same seed should produce identical null_hit_counts/null_tsd_lengths
    across two calls."""
    tsd_left, tsd_right = _make_flank_pair(flank_len=20, tsd_seq="GCATGCA")
    no_tsd_left, no_tsd_right = "A" * 20, "C" * 20
    left_seqs = [tsd_left] * 5 + [no_tsd_left] * 5
    right_seqs = [tsd_right] * 5 + [no_tsd_right] * 5

    hits1, lengths1 = shuffled_null_tsd_lengths(left_seqs, right_seqs, min_tsd=3, n_shuffles=15, seed=7)
    hits2, lengths2 = shuffled_null_tsd_lengths(left_seqs, right_seqs, min_tsd=3, n_shuffles=15, seed=7)

    assert hits1 == hits2
    assert lengths1 == lengths2


def test_tsd_significance_p_value_bounds():
    """p_value should always be in [0, 1]; enrichment should be inf when
    null_mean is 0, otherwise observed_hits / null_mean."""
    # left/right guarantee the shuffled-pairing null never finds a TSD
    # (an "A"-only flank can never share a run with a "C"-only flank), so
    # null_mean is forced to exactly 0 regardless of the seed.
    left_seqs = ["A" * 20] * 5
    right_seqs = ["C" * 20] * 5
    # tsd_results supplied independently of left_seqs/right_seqs above --
    # tsd_significance only uses them to regenerate the null, not to
    # recompute the observed hits.
    tsd_results = [{"acc_name": f"acc{i}", "tsd_len": 7 if i < 3 else 0} for i in range(5)]

    stats = tsd_significance(tsd_results, left_seqs, right_seqs, min_tsd=5, n_shuffles=10, seed=1)

    assert stats["observed_hits"] == 3
    assert stats["null_mean"] == 0.0
    assert stats["enrichment"] == float("inf")
    assert 0.0 <= stats["p_value"] <= 1.0


def test_build_split_entries_deduplicates_each_group():
    """Repeated (left, right) pairs should only appear once in tsd_entries /
    no_tsd_entries; repeated (left, padded_right) once in trimmed_entries."""
    flank_len = 10
    left, right = "AAAAAAAAAA", "XXXXXXXXXX"
    tsd_results = [
        {"acc_name": "acc1", "pairs": [(left, right)], "tsd_len": 4, "tsd_seq": "AAAA"},
        {"acc_name": "acc2", "pairs": [(left, right)], "tsd_len": 4, "tsd_seq": "AAAA"},  # identical pair
    ]

    tsd_entries, no_tsd_entries, trimmed_entries = build_split_entries(tsd_results, flank_len)
    assert len(tsd_entries) == 1
    assert len(trimmed_entries) == 1
    assert no_tsd_entries == []

    no_tsd_results = [
        {"acc_name": "acc3", "pairs": [("BBBBBBBBBB", "YYYYYYYYYY")], "tsd_len": 0, "tsd_seq": ""},
        {"acc_name": "acc4", "pairs": [("BBBBBBBBBB", "YYYYYYYYYY")], "tsd_len": 0, "tsd_seq": ""},
    ]
    tsd_entries2, no_tsd_entries2, trimmed_entries2 = build_split_entries(no_tsd_results, flank_len)
    assert len(no_tsd_entries2) == 1
    assert tsd_entries2 == [] and trimmed_entries2 == []


def test_build_split_entries_trimmed_padding_preserves_length():
    """Trimmed right flank should be right[tsd_len:] + 'N'*tsd_len, and its
    total length should equal flank_len (verify no off-by-one on padding)."""
    flank_len = 10
    left, right = "AAAAAAAAAA", "0123456789"

    for tsd_len in (1, 4, 9):  # small, middle, and near-full-length TSD
        tsd_results = [
            {"acc_name": f"acc_{tsd_len}", "pairs": [(left, right)], "tsd_len": tsd_len, "tsd_seq": right[:tsd_len]},
        ]
        _, _, trimmed_entries = build_split_entries(tsd_results, flank_len)
        assert len(trimmed_entries) == 1
        padded_right = trimmed_entries[0]["right"]
        assert len(padded_right) == flank_len
        assert padded_right == right[tsd_len:] + "N" * tsd_len


def test_build_split_entries_skips_accessions_with_no_valid_pairs():
    """An accession where every pair has 'none' on one side should be
    skipped (not added to any of the three groups)."""
    tsd_results = [
        {"acc_name": "acc1", "pairs": [("none", "CCCC"), ("AAAA", "none")], "tsd_len": 0, "tsd_seq": ""},
    ]
    tsd_entries, no_tsd_entries, trimmed_entries = build_split_entries(tsd_results, flank_len=4)
    assert tsd_entries == []
    assert no_tsd_entries == []
    assert trimmed_entries == []


def test_write_tsd_summary_percentages_sum_correctly(tmp_path):
    """With-TSD and without-TSD percentages in the written summary should sum
    to 100%."""
    tsd_results = (
        [{"acc_name": f"hit{i}", "tsd_len": 7, "tsd_seq": "GCATGCA"} for i in range(3)]
        + [{"acc_name": f"miss{i}", "tsd_len": 0, "tsd_seq": ""} for i in range(7)]
    )

    write_tsd_summary(tsd_results, "PCTTEST", str(tmp_path))

    summary = (tmp_path / "files" / "PCTTEST_TSD_SUMMARY.txt").read_text()
    assert "With TSD:         3 (30.0%)" in summary
    assert "Without TSD:      7 (70.0%)" in summary


def test_write_tsd_summary_empty_input_raises_zero_division(tmp_path):
    """KNOWN ISSUE: write_tsd_summary divides by `total` unguarded, so an
    empty tsd_results list raises ZeroDivisionError instead of writing an
    empty/placeholder summary. Documenting the current behaviour rather than
    silently masking it."""
    with pytest.raises(ZeroDivisionError):
        write_tsd_summary([], "EMPTY", str(tmp_path))


def test_find_tsds_detects_planted_tsd_at_correct_rate():
    """Build fake flank data where a known fraction of accessions carry a
    planted TSD of a known length, and the rest carry none, then confirm
    find_tsds recovers exactly that fraction and length -- not some other
    length, and not a suspiciously higher/lower hit rate from an accidental
    match in the "no TSD" filler sequences."""
    flank_len = 20
    tsd_len = 7
    n_with_tsd = 12
    n_without_tsd = 8

    tsd_left, tsd_right = _make_flank_pair(flank_len, "GCATGCA")
    no_tsd_left, no_tsd_right = "A" * flank_len, "C" * flank_len

    left_seqs = [tsd_left] * n_with_tsd + [no_tsd_left] * n_without_tsd
    right_seqs = [tsd_right] * n_with_tsd + [no_tsd_right] * n_without_tsd
    acc_names = [f"acc{i}" for i in range(n_with_tsd + n_without_tsd)]

    results = find_tsds(left_seqs, right_seqs, acc_names, min_tsd=3)

    hits = [r for r in results if r["tsd_len"] > 0]
    total = len(results)
    hit_rate = len(hits) / total

    assert len(hits) == n_with_tsd
    assert hit_rate == pytest.approx(n_with_tsd / total)
    assert all(r["tsd_len"] == tsd_len for r in hits)
    assert all(r["tsd_seq"] == "GCATGCA" for r in hits)


def test_run_tsd_analysis_end_to_end_smoke(tmp_path):
    """Same planted-TSD setup as test_find_tsds_detects_planted_tsd_at_correct_rate,
    but through the full run_tsd_analysis pipeline: checks the returned stats
    match the planted rate/length, TSD.txt and no_TSD.txt are written under
    <outdir>/files with the right number of (deduplicated) entries, and the
    summary/significance plots land under <outdir>/graphs."""
    flank_len = 20
    tsd_len = 7
    n_with_tsd = 12
    n_without_tsd = 8

    tsd_left, tsd_right = _make_flank_pair(flank_len, "GCATGCA")
    no_tsd_left, no_tsd_right = "A" * flank_len, "C" * flank_len

    left_seqs = [tsd_left] * n_with_tsd + [no_tsd_left] * n_without_tsd
    right_seqs = [tsd_right] * n_with_tsd + [no_tsd_right] * n_without_tsd
    acc_names = [f"acc{i}" for i in range(n_with_tsd + n_without_tsd)]

    stats = run_tsd_analysis(
        left_seqs, right_seqs, acc_names, flank_len, "TESTTSD", str(tmp_path),
        min_tsd=3, n_shuffles=20, seed=1,
    )

    total = n_with_tsd + n_without_tsd
    assert stats["observed_hits"] == n_with_tsd
    assert stats["obs_lengths"] == [tsd_len] * n_with_tsd

    files_dir = tmp_path / "files"
    graphs_dir = tmp_path / "graphs"
    # All planted TSD accessions share one (left, right) pair, so they
    # collapse to a single deduplicated entry; same for the no-TSD group.
    tsd_lines = (files_dir / "TESTTSD_TSD.txt").read_text().strip().splitlines()
    no_tsd_lines = (files_dir / "TESTTSD_no_TSD.txt").read_text().strip().splitlines()
    assert len(tsd_lines) == 1
    assert len(no_tsd_lines) == 1

    assert (files_dir / "TESTTSD_TSD_SUMMARY.txt").exists()
    assert (graphs_dir / "TESTTSD_TSD_NULL_DIST.png").exists()
    assert (graphs_dir / "TESTTSD_TSD_LENGTH_DIST.png").exists()


def test_run_tsd_analysis_seed_reproducibility(tmp_path):
    """Same seed should give identical null_hit_counts / p_value across two
    runs on the same input."""
    flank_len = 20
    tsd_left, tsd_right = _make_flank_pair(flank_len, "GCATGCA")
    no_tsd_left, no_tsd_right = "A" * flank_len, "C" * flank_len

    left_seqs = [tsd_left] * 6 + [no_tsd_left] * 4
    right_seqs = [tsd_right] * 6 + [no_tsd_right] * 4
    acc_names = [f"acc{i}" for i in range(10)]

    stats1 = run_tsd_analysis(
        left_seqs, right_seqs, acc_names, flank_len, "SEEDTEST", str(tmp_path / "run1"),
        min_tsd=3, n_shuffles=20, seed=99,
    )
    stats2 = run_tsd_analysis(
        left_seqs, right_seqs, acc_names, flank_len, "SEEDTEST", str(tmp_path / "run2"),
        min_tsd=3, n_shuffles=20, seed=99,
    )

    assert stats1["null_hit_counts"] == stats2["null_hit_counts"]
    assert stats1["p_value"] == stats2["p_value"]
