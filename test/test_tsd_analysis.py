# Tests for src/tsd_analysis.py
# (see also root-level tsd_test.py, which is a real-data sanity script, not a unit test)


def test_find_longest_tsd_exact_match():
    """A right flank starting with the same k bases that end the left flank
    should return the correct (longest) TSD length."""


def test_find_longest_tsd_below_min_tsd_returns_zero():
    """A shared suffix/prefix shorter than min_tsd should return 0, not the
    short match."""


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


def test_run_tsd_analysis_end_to_end_smoke(tmp_path):
    """Small synthetic left/right/acc_names input should produce TSD.txt and
    no_TSD.txt under <outdir>/files (trimmed entries are computed in memory
    by build_split_entries but no longer written to their own file), plus
    the TSD summary and significance plots under <outdir>/graphs, without
    raising."""


def test_run_tsd_analysis_seed_reproducibility(tmp_path):
    """Same seed should give identical null_hit_counts / p_value across two
    runs on the same input."""
