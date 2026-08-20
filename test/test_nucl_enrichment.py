# Tests for src/nucl_enrichment.py
#
# Only dinucleotide enrichment remains (trinucleotide enrichment was removed).
# The inner/outer split length is still hardcoded (7bp) rather than derived
# from TSD length + buffer -- worth testing current behaviour before changing it.


def test_count_dinucleotides_counts_all_overlapping_pairs():
    """A known sequence's dinucleotide counts should match manual counting,
    including overlaps (e.g. 'AAA' -> {'AA': 2})."""


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
    and <outdir>/graphs/<mge_name>_{LEFT,RIGHT}_DINUC_ENRICH.png."""
