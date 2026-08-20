# Tests for src/gc_content.py
#
# README says GC%/species analysis is being pulled out of the main pipeline,
# so this module's future is uncertain -- but while it exists it should be
# tested, especially the suspected bug below.


def test_gc_percent_known_sequence():
    """gc_percent('GCGC') == 100, gc_percent('ATAT') == 0,
    gc_percent('') == 0 (guard against division by zero)."""


def test_compute_gc_content_slices_do_not_depend_on_actual_flank_length():
    """SUSPECTED BUG: compute_gc_content hardcodes seq[:15] for left and
    seq[15:] for right regardless of the real flank length (no flank_len
    parameter exists on this function, unlike every other module in src/).
    Write a test with flanks longer than 15bp and confirm whether the GC%
    computed matches the whole flank or just an arbitrary 15bp window --
    this will show whether the hardcoding is a real bug before it's fixed."""


def test_compute_gc_content_deduplicates_sequences():
    """Duplicate left/right sequences should only be counted once each
    (list(set(...)) dedup) -- verify with a repeated sequence in the input."""


def test_compute_gc_content_summary_stats_match_input(tmp_path):
    """Written summary TSV should have correct mean/median/min/max/count per
    side, plus a 'total' row averaging left+right combined GC% (not the
    mean of the two per-side means, unless that's actually the same thing
    for equal-sized groups -- check the maths)."""


def test_compute_gc_content_writes_plot_and_tsv(tmp_path):
    """Should write <outdir>/graphs/gc_content_distribution.png and
    <outdir>/files/gc_content_summary.tsv without raising."""
