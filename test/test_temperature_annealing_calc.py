# Tests for src/temperature_annealing_calc.py
#
# Wired into flankalyser.py via run_melting_analysis() (--analyse_melting),
# and also usable standalone via main(). The Peyrard-Bishop breathing
# analysis was removed; only H-bond count and nearest-neighbor ΔG remain.
# Its shuffled control now goes through the shared src.utils.dinucleotide_shuffle
# (a proper dinucleotide-preserving shuffle) instead of a local mono-nucleotide one.


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


def test_run_melting_analysis_filters_n_containing_sequences():
    """Sequences containing any 'N' should be dropped before profiling
    (stricter than parse_input, which only drops all-N segments) -- confirm
    the logged usable-count reflects this."""


def test_run_melting_analysis_shuffled_control_is_dinucleotide_preserving():
    """The 'Shuffled' control line should come from src.utils.dinucleotide_shuffle
    -- verify the shuffled sequences have matching dinucleotide frequency
    counts to the originals (not just matching single-base composition)."""


def test_run_melting_analysis_end_to_end_smoke(tmp_path):
    """Small synthetic equal-length left/right sequences should produce
    <outdir>/graphs/<mge_name>_MELTING_PROFILE.png and a result dict with
    mean_hbonds_left/mean_hbonds_right, without raising."""


def test_run_melting_analysis_no_usable_sequences_returns_none():
    """If every sequence on one side contains an N, the function should log
    a warning and return None rather than raising on an empty list."""
