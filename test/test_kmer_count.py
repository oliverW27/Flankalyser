# Tests for src/kmer_count.py


def test_get_unique_kmers_counts_distinct_kmers_across_sequences():
    """Unique k-mers should be pooled across all input sequences (a k-mer
    repeated within/between sequences is only counted once)."""


def test_get_unique_kmers_k_equals_sequence_length():
    """When k == len(seq), exactly one k-mer per sequence should be produced."""


def test_get_unique_kmers_k_greater_than_sequence_length_returns_empty():
    """range(len(seq) - k + 1) goes negative -> should yield no k-mers, not
    raise an IndexError."""


def test_plot_unique_kmer_counts_both_output_path_is_not_outdir_relative():
    """KNOWN ISSUE: this function hardcodes 'graphs/kmer_distribution.png'
    instead of using an outdir/output_dir parameter like every other plotting
    function in this codebase. Test should assert current (buggy) behaviour
    first, then be updated once the output path is fixed to accept outdir."""


def test_plot_unique_kmer_counts_both_smoke(tmp_path, monkeypatch):
    """Run with small synthetic left/right sequences and default
    combined_props; assert it completes without raising and produces a plot
    (chdir into tmp_path to avoid polluting the repo, given the hardcoded
    path above)."""
