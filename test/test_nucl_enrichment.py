# Tests for src/nucl_enrichment.py
#
# README flags this module as needing work: inner/outer split lengths are
# hardcoded (7bp for dinucleotide, 25bp for trinucleotide) rather than
# derived from TSD length + buffer, and output file naming is described as
# "not great" -- worth testing current behaviour before changing it.


def test_count_dinucleotides_counts_all_overlapping_pairs():
    """A known sequence's dinucleotide counts should match manual counting,
    including overlaps (e.g. 'AAA' -> {'AA': 2})."""


def test_count_trinucleotides_counts_all_overlapping_triples():
    """Same as above for trinucleotides."""


def test_analyse_dinucleotide_enrichment_weighted_inner_outer_split():
    """For flank='left', inner should be the LAST inner_bp bases (near
    insertion site) and outer the first (flank_len - inner_bp); for
    flank='right' it should be reversed. Verify with a hand-built sequence
    where inner/outer have distinguishable content."""


def test_analyse_dinucleotide_enrichment_weighted_pseudocount_no_div_by_zero():
    """A dinucleotide absent from both inner and outer counts should still
    produce a finite log2 ratio (pseudocount avoids 0/0 or log(0))."""


def test_analyse_dinucleotide_enrichment_weighted_writes_expected_files(tmp_path):
    """Should write <outdir>/files/{left,right}_dinucleotide_enrichment.tsv
    and matching PNGs under <outdir>/graphs/."""


def test_analyse_trinucleotide_enrichment_weighted_top_n_selection():
    """Output should contain exactly top_n rows, sorted by descending
    absolute log2_enrichment."""


def test_analyse_trinucleotide_enrichment_weighted_inner_bp_exceeds_flank_len():
    """inner_bp=25 with a flank_len < 25 would make `split` negative --
    confirm what currently happens (slice with negative index) and whether
    that's the intended behaviour or a latent bug worth fixing."""
