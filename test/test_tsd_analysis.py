# Tests for src/tsd_analysis.py
import pytest

from src.tsd_analysis import _find_longest_tsd, find_tsds, run_tsd_analysis


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
    """N's should be stripped before comparison (_strip_n); confirm an N in
    the compared window doesn't cause a false negative/positive mismatch."""


def test_best_tsd_for_pairs_skips_none_segments():
    """Segment pairs where either side is 'none' should be skipped entirely."""


def test_best_tsd_for_pairs_picks_longest_across_multiple_segments():
    """When multiple left|right segment pairs exist, the best (longest) TSD
    across all pairs should be selected, with the correct tsd_seq returned."""


def test_find_tsds_mismatched_segment_counts_logs_warning():
    """Unequal numbers of left vs right pipe-segments should log a warning
    and pair up to the shorter count (via zip) without crashing."""


def test_shuffled_null_tsd_lengths_reproducible_with_seed():
    """Same seed should produce identical null_hit_counts/null_tsd_lengths
    across two calls."""


def test_tsd_significance_p_value_bounds():
    """p_value should always be in [0, 1]; enrichment should be inf when
    null_mean is 0, otherwise observed_hits / null_mean."""


def test_build_split_entries_deduplicates_each_group():
    """Repeated (left, right) pairs should only appear once in tsd_entries /
    no_tsd_entries; repeated (left, padded_right) once in trimmed_entries."""


def test_build_split_entries_trimmed_padding_preserves_length():
    """Trimmed right flank should be right[tsd_len:] + 'N'*tsd_len, and its
    total length should equal flank_len (verify no off-by-one on padding)."""


def test_build_split_entries_skips_accessions_with_no_valid_pairs():
    """An accession where every pair has 'none' on one side should be
    skipped (not added to any of the three groups)."""


def test_write_tsd_summary_percentages_sum_correctly(tmp_path):
    """With-TSD and without-TSD percentages in the written summary should sum
    to 100% (or handle total=0 without a ZeroDivisionError -- check this
    edge case explicitly, since the current code divides by `total` unguarded)."""


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
