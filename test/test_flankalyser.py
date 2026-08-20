# Tests for flankalyser.py CLI/pipeline logic.
#
# NOTE: flankalyser.py currently fails to import (see src.palindrome_analysis /
# src.tsd_detection imports at the top of the file, which don't exist in src/).
# These tests can't be written/run until that's fixed.


def test_parse_summary_file_parses_valid_rows():
    """acc/left/right TSV rows should map to three parallel lists in order."""


def test_parse_summary_file_skips_malformed_lines():
    """Lines with fewer than 3 tab-separated columns should be skipped (logged), not crash."""


def test_parse_summary_file_missing_file_exits():
    """A nonexistent summary_file path should log an error and sys.exit(1)."""


def test_extract_unique_sequences_splits_pipe_segments():
    """Multi-segment entries ('a|b|c') should be split on '|' and each segment
    considered independently for uniqueness."""


def test_extract_unique_sequences_skips_none_segments():
    """Segments equal to the literal string 'none' should be dropped."""


def test_extract_unique_sequences_enforces_flank_len():
    """Segments whose length != flank_len should be skipped with a warning,
    not silently included or crash the length-based counters downstream."""


def test_extract_unique_sequences_exclude_n_flag():
    """With exclude_n=True, segments containing 'N'/'n' should be excluded and counted."""


def test_extract_unique_sequences_rm_likely_snps_merges_single_mismatch():
    """With rm_likely_snps=True, sequences differing by exactly one Hamming-distance
    SNP should be merged into the most frequent variant; verify the correct
    variant (highest count) is kept and others are dropped."""


def test_extract_unique_sequences_rm_likely_snps_does_not_merge_two_mismatches():
    """Sequences differing by >1 SNP should NOT be merged."""


def test_load_gc_content_found_and_not_found():
    """load_gc_content should return the correct decimal GC% when the mge_name
    is present in the TSV, and None (with a warning) when it isn't."""


def test_create_default_proportions_shape():
    """Default proportions should be 0.25 for A/C/G/T, with a 'nan' sentinel
    at the joining position, length 2*flank_len + 1 per nucleotide."""


def test_generate_test_data_output_shape():
    """Randomly generated test data should have num_accs rows, each with a
    left and right flank of length flank_len."""


def test_main_end_to_end_smoke():
    """Once imports are fixed: run main() (or an equivalent pipeline call) on
    the small example in data/flanks_ISFinder_IS1F.txt and assert output files/
    graphs are created under --outdir without raising."""
