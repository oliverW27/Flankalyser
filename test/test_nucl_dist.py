# Tests for src/nucl_dist.py


def test_count_nucl_counts_each_position_correctly():
    """A small hand-built list of equal-length sequences should produce exact
    per-position A/C/G/T counts for both left and right."""


def test_count_nucl_warns_on_length_mismatch():
    """A sequence whose length != FLANK_LEN should still be counted for the
    positions it has but log a warning (verify it doesn't crash/misalign)."""


def test_count_nucl_ignores_n_but_warns_on_other_unexpected_chars():
    """'N' should be silently skipped; any other non-ACGT character should log
    a warning and not increment any bucket."""


def test_nucl_dist_stat_sig_default_25pct_expectation():
    """With gc_content=None, expected proportions should be 25% each; verify
    p-values against a known chi-square case (e.g. perfectly uniform counts
    -> p close to 1, heavily skewed counts -> p close to 0)."""


def test_nucl_dist_stat_sig_gc_adjusted_expectation():
    """With gc_content=0.6, expected G/C should each be 0.3 and A/T each 0.2;
    verify this changes the resulting p-values vs. the default case."""


def test_nucl_dist_stat_sig_handles_zero_counts_at_a_position():
    """A position with zero total counts should not crash chisquare (currently
    caught by try/except -> NaN); confirm NaN is returned, not an exception."""


def test_calc_nucl_props_combines_left_and_right_with_nan_separator():
    """Output should be left proportions + [nan] + right proportions, each
    position's per-nucleotide proportions summing to ~1 (or NaN if no counts)."""


def test_calc_nucl_props_all_zero_counts_gives_nan():
    """A locus with zero counts across all four nucleotides should yield NaN
    proportions rather than a ZeroDivisionError."""


def test_adjust_gc_content_baseline_calculation():
    """Baseline should be the mean proportion in the outer `outer_bp` positions
    of each flank (excluding the NaN join position); verify adjusted values
    equal original - (baseline - 0.25)."""


def test_adjust_gc_content_missing_outer_data_defaults_to_025():
    """If the outer region has no valid (non-NaN) data for a nucleotide, the
    baseline should default to 0.25 with a warning, not crash."""


def test_get_stats_ast_bonferroni_thresholds():
    """Verify asterisk assignment respects Bonferroni-corrected thresholds
    (0.05/n, 0.01/n, 0.001/n) and that 'NaN' entries map to a blank space,
    not an exception from comparing str < float."""


def test_graph_smoke_writes_png(tmp_path):
    """graph() should write a PNG to <outdir>/graphs/DNA_freq_test_<name>.png
    without raising, for both gc_content=None and a supplied gc_content."""
