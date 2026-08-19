"""
Empirical comparison tests: palindrome v1 (short script) vs v2 (palindrome_analysis.py).

These tests don't just check "pass/fail" — they build small synthetic TSVs that
isolate one behavioral difference at a time, run both scripts, and assert on
what actually differs. Read the docstring on each test for what it's proving.

Fill in the globals below, then run:
    pytest test_palindrome_compare.py -v -s
"""
import importlib.util
import sys
import types
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# FILL THESE IN
# ---------------------------------------------------------------------------
V1_SCRIPT_PATH = "src/palindrome_analysis.py"              # the short script
V2_SCRIPT_PATH = "src/find_pals.py"   # the longer / attached script
REAL_FLANKS_TSV = "../flankalyser/IS3/flanks_ISFinder_IS1397.txt"                                # optional: real acc/left/right TSV
                                                     # for a side-by-side sanity check;
                                                     # leave "" to skip that test
MGE_NAME = "test_mge"
GRAPH_OUTPUT_DIR = "."
# ---------------------------------------------------------------------------


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_tsv(path, rows):
    """rows: list of (accession, left, right) tuples."""
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")


def _load_module(path, name):
    if not Path(path).exists():
        pytest.skip(f"{name} not found at {path!r} — set the path at the top of this file")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def v1():
    # v1 imports `from src.utils import generate_position_weighted_random_seq`.
    # If that package isn't importable in this environment, stub it so the
    # module still loads — the stub gets overridden per-test via monkeypatch.
    try:
        import src.utils  # noqa: F401
    except ImportError:
        stub_pkg = types.ModuleType("src")
        stub_utils = types.ModuleType("src.utils")
        stub_utils.generate_position_weighted_random_seq = (
            lambda combined_props, start, end: "A" * (end - start)
        )
        sys.modules["src"] = stub_pkg
        sys.modules["src.utils"] = stub_utils
    return _load_module(V1_SCRIPT_PATH, "palindrome_v1")


@pytest.fixture(scope="session")
def v2():
    return _load_module(V2_SCRIPT_PATH, "palindrome_v2")


@pytest.fixture
def run_v1(monkeypatch, v1):
    """
    Calls v1's process_summary_file_centered_palindromes and returns
    (results, random_results) by intercepting the call to plot_palindrome_results
    instead of letting it write a PNG.
    """
    def _run(tsv_path, min_len=5, max_len=8):
        captured = {}

        def fake_plot(results, random_results, output_path=None):
            captured["results"] = dict(results)
            captured["random_results"] = dict(random_results)

        monkeypatch.setattr(v1, "plot_palindrome_results", fake_plot)
        monkeypatch.setattr(
            v1, "generate_position_weighted_random_seq",
            lambda combined_props, start, end: "A" * (end - start),
            raising=False,
        )
        v1.process_summary_file_centered_palindromes(
            left_seqs=None, right_seqs=None, combined_props=None,
            path=str(tsv_path), min_len=min_len, max_len=max_len,
        )
        return captured["results"], captured["random_results"]

    return _run


# ---------------------------------------------------------------------------
# 1. Sanity check — a clean palindrome sitting exactly at the L/R junction,
#    with matching flank lengths (50bp), should be picked up by both.
# ---------------------------------------------------------------------------
def test_boundary_palindrome_detected_by_both(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"          # 6bp arm
    rc_arm = revcomp(arm)   # "GTACGT"
    left = "G" * 44 + arm        # 50bp, arm right at the end
    right = rc_arm + "G" * 44    # 50bp, arm right at the start
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    assert results[5] == 1 and results[6] == 1  # nested: 6bp arm also satisfies k=5
    assert results[7] == 0 and results[8] == 0

    seqs = v2.load_sequences(str(tsv))
    acc, seq = seqs[0]
    hit = v2.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=50
    )
    assert hit["arm_len"] == 6


# ---------------------------------------------------------------------------
# 2. Counting scheme difference: v1 buckets the SAME hit into every k it
#    satisfies (cumulative), v2 buckets it once, at its true longest length.
#    A single 6bp palindrome shows up as TWO hits in v1 (k=5 and k=6) but
#    ONE hit in v2 (arm_len=6).
# ---------------------------------------------------------------------------
def test_v1_counts_cumulatively_v2_counts_longest_only(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    v1_hit_count_across_ks = sum(results.values())
    assert v1_hit_count_across_ks == 2  # counted at both k=5 and k=6

    seqs = v2.load_sequences(str(tsv))
    acc, seq = seqs[0]
    hit = v2.find_longest_inverted_repeat(seq, max_spacer=0, min_arm=5, max_mismatches=0)
    # v2 reports one record for this sequence, with the true arm length
    assert hit["arm_len"] == 6


# ---------------------------------------------------------------------------
# 3. Multi-segment ("|") rows: v1 iterates over every left|right pair, v2
#    skips the row entirely.
# ---------------------------------------------------------------------------
def test_pipe_separated_multi_segment_rows(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    # two segments per side, pipe-separated
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", f"{left}|{left}", f"{right}|{right}")])

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    # v1 checks all left x right combinations -> 2x2 = 4 pair comparisons,
    # each satisfying both k=5 and k=6
    assert results[5] == 4 and results[6] == 4

    seqs = v2.load_sequences(str(tsv))
    assert seqs == []  # v2's load_sequences drops any row containing "|"


# ---------------------------------------------------------------------------
# 4. Search scope: v1 only ever looks exactly at the L/R boundary. v2's
#    default (non-centered) mode searches the whole concatenated sequence,
#    so it finds inverted repeats that sit well away from the junction —
#    v1 has no way to see these at all.
# ---------------------------------------------------------------------------
def test_off_junction_palindrome_only_found_by_v2(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 20                              # no signal here
    right = "G" * 7 + arm + rc_arm + "G" * 7      # IR sits inside the right flank,
                                                   # nowhere near the boundary
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    assert sum(results.values()) == 0  # boundary itself (G|G) has no complementarity

    seqs = v2.load_sequences(str(tsv))
    acc, seq = seqs[0]
    hit = v2.find_longest_inverted_repeat(seq, max_spacer=0, min_arm=3, max_mismatches=0)
    assert hit["arm_len"] == 6  # found it, away from the junction


# ---------------------------------------------------------------------------
# 5. v2's --centered mode hardcodes the palindrome centre at index 50. That's
#    only correct if each flank is actually 50bp. With shorter flanks it
#    silently finds nothing, rather than erroring — worth knowing before
#    trusting --centered results on data with a different flank length.
# ---------------------------------------------------------------------------
def test_v2_centered_mode_assumes_50bp_flanks(tmp_path, v2):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 20                              # 26bp right flank -> seq len 46, not 100
    right = "G" * 7 + arm + rc_arm + "G" * 7
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    seqs = v2.load_sequences(str(tsv))
    acc, seq = seqs[0]
    assert len(seq) == 46

    hit = v2.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=50
    )
    assert hit["arm_len"] == 0  # centre=50 is out of range for a 46bp sequence


# ---------------------------------------------------------------------------
# 6. Mismatch tolerance: v2 can allow mismatches within an arm (max_mismatches),
#    v1 has no such concept — any imperfection anywhere in the compared window
#    fails that k entirely, it doesn't "count mismatches" at all.
# ---------------------------------------------------------------------------
def test_mismatch_tolerance_v2_only(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"
    rc_arm_perfect = revcomp(arm)          # "GTACGT"
    rc_arm_1_mismatch = "GTACGA"           # last base wrong vs. perfect revcomp
    left = "G" * 44 + arm
    right = rc_arm_1_mismatch + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right)])

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    assert results[6] == 0  # exact 6bp comparison fails because of the mismatch

    seqs = v2.load_sequences(str(tsv))
    acc, seq = seqs[0]
    exact = v2.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=0, fixed_centre=50
    )
    tolerant = v2.find_longest_inverted_repeat(
        seq, max_spacer=0, min_arm=3, max_mismatches=1, fixed_centre=50
    )
    assert exact["arm_len"] < 6         # same shortfall v1 sees
    assert tolerant["arm_len"] == 6     # v2 recovers the full arm by allowing 1 mismatch


# ---------------------------------------------------------------------------
# 7. Deduplication: v2's load_sequences drops duplicate concatenated
#    sequences, v1 has no dedup — repeated rows just multiply the count.
# ---------------------------------------------------------------------------
def test_duplicate_sequences_deduped_only_by_v2(tmp_path, v1, v2, run_v1):
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 44 + arm
    right = rc_arm + "G" * 44
    tsv = tmp_path / "flanks.tsv"
    write_tsv(tsv, [("acc1", left, right), ("acc2", left, right)])  # identical seqs

    results, _ = run_v1(tsv, min_len=5, max_len=8)
    assert results[6] == 2  # counted once per row

    seqs = v2.load_sequences(str(tsv))
    assert len(seqs) == 1  # deduped down to one unique sequence


# ---------------------------------------------------------------------------
# 8. Optional: run both on the same real dataset and print a side-by-side
#    summary. No strict assertions here since real data outcomes aren't
#    known in advance — this is for eyeballing overall differences.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not REAL_FLANKS_TSV, reason="set REAL_FLANKS_TSV to enable")
def test_real_data_side_by_side(v1, v2, run_v1):
    results, _ = run_v1(REAL_FLANKS_TSV, min_len=5, max_len=8)
 
    seqs = v2.load_sequences(REAL_FLANKS_TSV)
    obs_arm_lens, null_arm_lens = [], []
    null_hits = []
    min_arm, max_spacer, max_mismatches, n_shuffles = 10, 1, 0, 30
 
    for _, seq in seqs:
        hit = v2.find_longest_inverted_repeat(
            seq, max_spacer=max_spacer, min_arm=min_arm, max_mismatches=max_mismatches
        )
        if hit["arm_len"] >= min_arm:
            obs_arm_lens.append(hit["arm_len"])
    observed_hits = len(obs_arm_lens)
 
    for _ in range(n_shuffles):
        count = 0
        for _, seq in seqs:
            shuffled = v2.dinucleotide_shuffle(seq)
            hit = v2.find_longest_inverted_repeat(
                shuffled, max_spacer=max_spacer, min_arm=min_arm, max_mismatches=max_mismatches
            )
            if hit["arm_len"] >= min_arm:
                count += 1
                null_arm_lens.append(hit["arm_len"])
        null_hits.append(count)
 
    p_val = sum(n >= observed_hits for n in null_hits) / n_shuffles
 
    print("\n--- v1 (per-k, cumulative, boundary-only) ---")
    for k, count in sorted(results.items()):
        print(f"  k={k}: {count} hits")
    print(f"\n--- v2 (min_arm={min_arm}, boundary or anywhere, spacer<={max_spacer}) ---")
    print(f"  {observed_hits}/{len(seqs)} sequences with an arm >= {min_arm}bp")
 
    out_dir = Path(GRAPH_OUTPUT_DIR)
    out_dir.mkdir(parents=True, exist_ok=True)
 
    null_dist_path = out_dir / f"{MGE_NAME}_null_distribution.png"
    v2.plot_null_distribution(null_hits, observed_hits, p_val, MGE_NAME, null_dist_path)
    print(f"\nNull distribution plot: {null_dist_path.resolve()}")
 
    arm_dist_path = out_dir / f"{MGE_NAME}_arm_length_distribution.png"
    v2.plot_arm_length_distribution(obs_arm_lens, null_arm_lens, n_shuffles, MGE_NAME, arm_dist_path)
    print(f"Arm length distribution plot: {arm_dist_path.resolve()}")
