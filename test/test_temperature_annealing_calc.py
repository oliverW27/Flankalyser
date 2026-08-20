# Tests for src/temperature_annealing_calc.py
# Standalone script, not currently imported by flankalyser.py or anything else
# in src/ -- confirm with the user whether it's meant to be wired into the
# main pipeline or kept as a separate tool before investing heavily here.


def test_parse_input_splits_pipe_segments_and_dedupes():
    """Left/right pipe-separated segments should become independent entries
    in the lefts/rights sets, deduplicated, uppercased."""


def test_parse_input_excludes_all_n_sequences():
    """A segment that is entirely 'N' characters should be excluded; a
    segment with SOME Ns should still be included (only all-N is filtered
    here, unlike find_pals.py which rejects any internal N)."""


def test_hbond_per_seq_known_sequence():
    """hbond_per_seq('ATGC') == [2, 2, 3, 3]; unknown characters should map
    to NaN, not raise a KeyError."""


def test_dg_per_seq_known_dinucleotides():
    """dg_per_seq should return NN_DG values for each overlapping dinucleotide
    step; an unknown dinucleotide should produce NaN."""


def test_profile_with_sem_handles_variable_length_sequences():
    """Sequences of different lengths should be right-padded with NaN before
    nanmean/nanstd, not truncated or raise a shape error."""


def test_profile_with_sem_empty_input_returns_none():
    """Empty seqs list (or all per_seq_fn results None) should return
    (None, None) rather than raising."""


def test_pb_breathing_profile_uses_most_common_length_only():
    """Sequences whose length doesn't match the modal length in the input
    should be silently excluded before the PB calculation -- confirm this
    with a mixed-length input and check the returned profile length matches
    the modal length, not the max or min."""


def test_pb_breathing_profile_max_seqs_subsampling_reproducible():
    """With max_seqs set and a fixed internal seed (42), repeated calls on
    the same input should subsample the same sequences."""


def test_shuffle_seqs_preserves_composition_changes_order():
    """Shuffled sequence should have identical character multiset to input,
    but (with high probability, for a long-enough sequence) a different
    order."""
