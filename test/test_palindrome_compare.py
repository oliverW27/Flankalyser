"""
Behavioral tests for src/find_pals.py's palindrome / inverted-repeat detection.

Originally compared find_pals.py ("v2") against an older inline palindrome
script ("v1", src/palindrome_analysis.py) that no longer exists in this repo
-- find_pals.py is the sole, superseding implementation now, so these tests
exercise it directly instead of diffing it against something that's gone.
Read the docstring on each test for the behaviour it's pinning down.

Run:
    pytest test_palindrome_compare.py -v -s
"""
from pathlib import Path

import pytest

from src import find_pals

REAL_FLANKS_TSV = "data/flanks_ISFinder_IS1F.txt"
MGE_NAME = "test_mge"


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_tsv(path, rows):
    """rows: list of (accession, left, right) tuples."""
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# 1. A clean palindrome sitting exactly at the L/R junction, in centered mode.
# ---------------------------------------------------------------------------
def test_boundary_palindrome_detected_with_correct_arm_length(tmp_path):
    arm = "ACGTAC"              # 6bp arm
    rc_arm = revcomp(arm)       # "GTACGT"
    left = "G" * 44 + arm       # 50bp, arm right at the end
    right = rc_arm + "G" * 44   # 50bp, arm right at the start
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = find_pals.load_sequences(str(tsv))
    acc, seq = seqs[0]
    assert len(seq) == 100

    hit = find_pals.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=len(seq) // 2
    )
    assert hit["arm_len"] == 6


# ---------------------------------------------------------------------------
# 2. find_longest_inverted_repeat reports the longest arm for a sequence,
#    not one hit per length threshold it happens to satisfy.
# ---------------------------------------------------------------------------
def test_reports_longest_arm_not_every_satisfying_length(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = find_pals.load_sequences(str(tsv))
    acc, seq = seqs[0]
    hit = find_pals.find_longest_inverted_repeat(seq, max_spacer=0, min_arm=5, max_mismatches=0)
    assert hit["arm_len"] == 6


# ---------------------------------------------------------------------------
# 3. Multi-segment ("|") rows are dropped entirely, not expanded.
# ---------------------------------------------------------------------------
def test_multi_segment_rows_are_skipped_entirely(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", f"{left}|{left}", f"{right}|{right}")])

    seqs = find_pals.load_sequences(str(tsv))
    assert seqs == []  # any row containing "|" is dropped


# ---------------------------------------------------------------------------
# 4. Default (non-centered) mode searches the whole concatenated sequence,
#    so it finds inverted repeats that sit well away from the L/R junction.
# ---------------------------------------------------------------------------
def test_off_junction_palindrome_found_when_not_centered(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 20                              # no signal here
    right = "G" * 7 + arm + rc_arm + "G" * 7      # IR sits inside the right flank,
                                                   # nowhere near the boundary
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = find_pals.load_sequences(str(tsv))
    acc, seq = seqs[0]
    hit = find_pals.find_longest_inverted_repeat(seq, max_spacer=0, min_arm=3, max_mismatches=0)
    assert hit["arm_len"] == 6


# ---------------------------------------------------------------------------
# 5. fixed_centre restricts the search to that exact index -- a palindrome
#    centred elsewhere is invisible to it, even though the same sequence is
#    found once the search is unrestricted. (This is why
#    run_palindrome_analysis derives fixed_centre from the sequence's own
#    midpoint instead of a hardcoded value -- see flankalyser_test.py's
#    regression test for that fix.)
# ---------------------------------------------------------------------------
def test_centered_mode_only_tests_the_exact_fixed_index(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 20
    right = "G" * 7 + arm + rc_arm + "G" * 7   # IR centred inside the right flank
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = find_pals.load_sequences(str(tsv))
    acc, seq = seqs[0]

    boundary_only = find_pals.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=len(left)
    )
    assert boundary_only["arm_len"] == 0  # nothing sits at the L/R boundary itself

    unrestricted = find_pals.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0
    )
    assert unrestricted["arm_len"] == 6  # found once the whole sequence is searched


# ---------------------------------------------------------------------------
# 6. Mismatch tolerance: an arm with one wrong base can still be recovered
#    by allowing max_mismatches=1; max_mismatches=0 stops short at the flaw.
# ---------------------------------------------------------------------------
def test_mismatch_tolerance_recovers_the_full_arm(tmp_path):
    arm = "ACGTAC"
    rc_arm_1_mismatch = "GTACGA"   # last base wrong vs. the perfect revcomp "GTACGT"
    left = "G" * 44 + arm
    right = rc_arm_1_mismatch + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = find_pals.load_sequences(str(tsv))
    acc, seq = seqs[0]
    centre = len(seq) // 2

    exact = find_pals.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=centre
    )
    tolerant = find_pals.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=1, fixed_centre=centre
    )
    assert exact["arm_len"] < 6
    assert tolerant["arm_len"] == 6


# ---------------------------------------------------------------------------
# 7. load_sequences deduplicates identical concatenated sequences.
# ---------------------------------------------------------------------------
def test_duplicate_sequences_deduplicated(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right), ("acc2", left, right)])  # identical seqs

    seqs = find_pals.load_sequences(str(tsv))
    assert len(seqs) == 1


# ---------------------------------------------------------------------------
# 8. The shared dinucleotide shuffle (src.utils, re-exported here) should
#    preserve length and single-base composition, only reorder the sequence.
# ---------------------------------------------------------------------------
def test_dinucleotide_shuffle_preserves_composition_and_length():
    seq = "ACGTACGTACGTACGTGGGGCCCCAAAATTTT"
    shuffled = find_pals.dinucleotide_shuffle(seq)
    assert len(shuffled) == len(seq)
    assert sorted(shuffled) == sorted(seq)


# ---------------------------------------------------------------------------
# 9. End-to-end smoke test: a planted palindrome should show up as an
#    observed hit, and the summary/plot files should all be written.
# ---------------------------------------------------------------------------
def test_run_palindrome_analysis_end_to_end_smoke(tmp_path):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    # A fixed pseudo-random 50/50bp pair with no long inverted repeat of its
    # own (checked: max arm_len is 3 at min_arm=1), so it's a clean "no hit"
    # control alongside the planted palindrome above.
    control_left = "AACATAGCTACCAACAGACCTAATATGCTTCTAGGTCCTGGGGTCGCATA"
    control_right = "TTGACCACGCAAGACAGGTGTTGACATCTAAACCCCCAACACAGAGCTCC"
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right), ("acc2", control_left, control_right)])

    find_pals.run_palindrome_analysis(
        input_tsv=str(tsv), mge_name=MGE_NAME, output_dir=str(tmp_path),
        min_arm=6, max_spacer=0, max_mismatches=0, n_shuffles=20, seed=1,
    )

    summary_path = tmp_path / "files" / f"{MGE_NAME}_PALINDROME_SUMMARY.txt"
    assert summary_path.exists()
    assert "Observed hits:        1/2" in summary_path.read_text()
    assert (tmp_path / "graphs" / f"{MGE_NAME}_PALINDROME_NULL_DIST.png").exists()
    assert (tmp_path / "graphs" / f"{MGE_NAME}_PALINDROME_ARM_LEN_DIST.png").exists()


# ---------------------------------------------------------------------------
# 10. Optional: run on the real example dataset and print the summary for a
#     manual eyeball. No strict assertions since real-data outcomes aren't
#     known in advance.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not Path(REAL_FLANKS_TSV).exists(), reason=f"{REAL_FLANKS_TSV!r} not found")
def test_real_data_smoke(tmp_path):
    find_pals.run_palindrome_analysis(
        input_tsv=REAL_FLANKS_TSV, mge_name="IS1F_REAL", output_dir=str(tmp_path),
        min_arm=6, max_spacer=2, max_mismatches=0, n_shuffles=20, seed=1,
    )
    summary_path = tmp_path / "files" / "IS1F_REAL_PALINDROME_SUMMARY.txt"
    assert summary_path.exists()
    print("\n" + summary_path.read_text())
