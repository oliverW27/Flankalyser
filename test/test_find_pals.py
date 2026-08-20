# Tests for src/find_pals.py
# (see also root-level palindrome_test.py, which compares this module against a
# now-missing src/palindrome_analysis.py "v1" -- those tests currently skip)


def test_find_longest_inverted_repeat_perfect_palindrome():
    """A sequence with a known perfect inverted repeat at a known centre should
    return the correct arm_len, spacer_size, and arm_seq."""


def test_find_longest_inverted_repeat_respects_min_arm():
    """An inverted repeat shorter than min_arm should not be reported (arm_len
    stays 0 / below best threshold)."""


def test_find_longest_inverted_repeat_mismatch_tolerance():
    """With max_mismatches=1, an arm with exactly one mismatch should be found
    at full length; with max_mismatches=0 it should stop short at the
    mismatch position."""


def test_find_longest_inverted_repeat_spacer_range():
    """Palindromes with spacers of varying size (0..max_spacer) should all be
    detected; a spacer larger than max_spacer should not be."""


def test_find_longest_inverted_repeat_fixed_centre_restricts_search():
    """With fixed_centre set, only that index should be tried -- a palindrome
    centred elsewhere in the sequence should NOT be found even though it
    would be found in unrestricted (fixed_centre=None) mode."""


def test_find_longest_inverted_repeat_search_window_bounds():
    """search_start/search_end should restrict where arms can extend to;
    an inverted repeat that would extend past the window boundary should be
    truncated or excluded, not read out-of-bounds."""


def test_dinucleotide_shuffle_preserves_dinucleotide_frequencies():
    """Shuffled sequence should have identical dinucleotide frequency counts
    to the input, only order changed."""


def test_dinucleotide_shuffle_preserves_length_and_composition():
    """Output length and single-nucleotide composition should match input."""


def test_dinucleotide_shuffle_raises_after_max_attempts_on_unshuffleable_seq():
    """A pathological sequence where the shuffle can't complete an Eulerian
    path should raise RuntimeError after max_attempts, not hang or silently
    return a short sequence."""


def test_load_sequences_skips_multisegment_rows():
    """Rows where left or right contains '|' should be dropped entirely."""


def test_load_sequences_strips_trailing_ns_from_right_and_leading_from_left():
    """N-padding from TSD trimming should be stripped before concatenation;
    verify the exact strip logic (right.rstrip, left.rstrip+lstrip)."""


def test_load_sequences_rejects_internal_ns():
    """A sequence with an N in the middle (post-strip) should be excluded."""


def test_load_sequences_deduplicates_identical_concatenated_sequences():
    """Two accessions producing the same concatenated left+right sequence
    should only appear once in the output."""


def test_load_sequences_respects_min_len():
    """Concatenated sequences shorter than min_len should be excluded."""


def test_run_palindrome_analysis_end_to_end_smoke(tmp_path):
    """Small synthetic TSV with a known palindrome should produce a summary
    txt and two PNGs under outdir, with observed_hits reflecting the known
    planted palindrome."""


def test_run_palindrome_analysis_seed_reproducibility(tmp_path):
    """Same seed should give identical p_val / null_hits across two runs."""
