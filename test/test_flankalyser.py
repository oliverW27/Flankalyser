# Tests for flankalyser.py CLI/pipeline logic.


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


def test_main_end_to_end_smoke():
    """Run main() (or an equivalent pipeline call) on the small example in
    data/flanks_ISFinder_IS1F.txt and assert output files/graphs are created
    under --outdir without raising."""


def test_main_palindrome_centered_uses_flank_boundary_not_hardcoded_50():
    """--analyse_palindromes --palindrome_centered should use the actual
    left/right flank boundary (total_len // 2) as the palindrome centre,
    not a hardcoded index -- run with flank_len != 50 and confirm it still
    finds a planted boundary palindrome."""
