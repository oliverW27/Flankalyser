# Tests for src/logo.py


def test_seqs_to_counts_builds_correct_position_frequency_matrix():
    """Given a small hand-built list of equal-length sequences, the resulting
    count matrix (rows=bases, cols=positions before transpose) should match
    manual counts exactly."""


def test_generate_sequence_logo_filters_variable_length_sequences():
    """Sequences whose length differs from the first sequence in the list
    should be silently dropped (per the current filtering logic) -- confirm
    this is intended, since it could silently discard most of the data if
    the first sequence happens to be an outlier length."""


def test_generate_sequence_logo_empty_input_after_filtering():
    """If left_seqs or right_seqs is empty (or becomes empty after length
    filtering), confirm what currently happens -- likely an IndexError on
    left_seqs[0] before filtering even runs. Worth an explicit test to
    document the failure mode."""


def test_generate_sequence_logo_writes_expected_png_files(tmp_path):
    """Should write <outdir>/graphs/left_flank_sequence_logo.png and
    right_flank_sequence_logo.png without raising, for valid equal-length
    input."""
