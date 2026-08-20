"""
Tests for flankalyser.py -- the main CLI/pipeline entry point.

Runnable standalone from the repo root:
    pytest flankalyser_test.py -v
"""
import logging
import sys

import pytest

import flankalyser as fk


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_tsv(path, rows):
    """rows: list of (accession, left, right) tuples."""
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(row) + "\n")


@pytest.fixture
def logger():
    return logging.getLogger("flankalyser_test")


def run_main(monkeypatch, args):
    """Invoke flankalyser.main() as if called from the CLI with `args`."""
    monkeypatch.setattr(sys, "argv", ["flankalyser.py"] + args)
    fk.main()


# ---------------------------------------------------------------------------
# parse_summary_file
# ---------------------------------------------------------------------------

def test_parse_summary_file_parses_valid_rows(tmp_path, logger):
    tsv = tmp_path / "summary.txt"
    write_tsv(tsv, [("acc1", "ACGT", "TGCA"), ("acc2", "GGGG", "CCCC")])

    l_data, r_data, acc_names = fk.parse_summary_file(str(tsv), logger)

    assert acc_names == ["acc1", "acc2"]
    assert l_data == ["ACGT", "GGGG"]
    assert r_data == ["TGCA", "CCCC"]


def test_parse_summary_file_skips_malformed_lines(tmp_path, logger):
    tsv = tmp_path / "summary.txt"
    tsv.write_text("acc1\tACGT\tTGCA\nnot_enough_columns\nacc2\tGGGG\tCCCC\n")

    l_data, r_data, acc_names = fk.parse_summary_file(str(tsv), logger)

    assert acc_names == ["acc1", "acc2"]
    assert len(l_data) == len(r_data) == 2


def test_parse_summary_file_missing_file_exits(tmp_path, logger):
    missing = tmp_path / "does_not_exist.txt"
    with pytest.raises(SystemExit):
        fk.parse_summary_file(str(missing), logger)


# ---------------------------------------------------------------------------
# extract_unique_sequences
# ---------------------------------------------------------------------------

def test_extract_unique_sequences_splits_pipe_segments(logger):
    l_data = ["ACGTACGTAC|TGCATGCATG"]
    r_data = ["GGGGGGGGGG|CCCCCCCCCC"]
    acc_names = ["acc1"]

    l_uniq, r_uniq, l_names, r_names = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger
    )

    assert set(l_uniq) == {"ACGTACGTAC", "TGCATGCATG"}
    assert set(r_uniq) == {"GGGGGGGGGG", "CCCCCCCCCC"}


def test_extract_unique_sequences_skips_none_segments(logger):
    l_data = ["ACGTACGTAC|none"]
    r_data = ["none|CCCCCCCCCC"]
    acc_names = ["acc1"]

    l_uniq, r_uniq, l_names, r_names = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger
    )

    assert l_uniq == ["ACGTACGTAC"]
    assert r_uniq == ["CCCCCCCCCC"]


def test_extract_unique_sequences_enforces_flank_len(logger):
    l_data = ["ACGT", "ACGTACGTAC"]  # first is the wrong length for flank_len=10
    r_data = ["TGCA", "TGCATGCATG"]
    acc_names = ["short_acc", "good_acc"]

    l_uniq, r_uniq, l_names, r_names = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger
    )

    assert l_uniq == ["ACGTACGTAC"]
    assert l_names == ["good_acc"]


def test_extract_unique_sequences_exclude_n_flag(logger):
    l_data = ["ACGTACGTAN", "ACGTACGTAC"]
    r_data = ["TGCATGCATG", "TGCATGCATG"]
    acc_names = ["has_n", "clean"]

    included_l, _, _, _ = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger, exclude_n=False
    )
    assert "ACGTACGTAN" in included_l

    excluded_l, _, _, _ = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger, exclude_n=True
    )
    assert "ACGTACGTAN" not in excluded_l
    assert "ACGTACGTAC" in excluded_l


def test_extract_unique_sequences_deduplicates(logger):
    l_data = ["ACGTACGTAC", "ACGTACGTAC"]  # identical sequence, two accessions
    r_data = ["TGCATGCATG", "TGCATGCATG"]
    acc_names = ["acc1", "acc2"]

    l_uniq, r_uniq, l_names, r_names = fk.extract_unique_sequences(
        l_data, r_data, acc_names, flank_len=10, logger=logger
    )

    assert l_uniq == ["ACGTACGTAC"]
    assert r_uniq == ["TGCATGCATG"]


# ---------------------------------------------------------------------------
# main() -- end to end, against the real example data
# ---------------------------------------------------------------------------

def test_main_default_run_produces_composition_and_enrichment_outputs(tmp_path, monkeypatch):
    run_main(monkeypatch, [
        "--mge_name", "IS1F",
        "--summary_file", "data/flanks_ISFinder_IS1F.txt",
        "--flank_len", "50",
        "--outdir", str(tmp_path),
    ])

    graphs = {p.name for p in (tmp_path / "graphs").glob("*.png")}
    assert "IS1F_DNA_FREQ.png" in graphs
    assert "IS1F_LEFT_LOGO.png" in graphs
    assert "IS1F_RIGHT_LOGO.png" in graphs
    assert "IS1F_LEFT_DINUC_ENRICH.png" in graphs
    assert "IS1F_RIGHT_DINUC_ENRICH.png" in graphs

    # TSD/palindrome/bendability/melting weren't requested -- their outputs
    # shouldn't exist.
    assert "IS1F_TSD_NULL_DIST.png" not in graphs
    assert "IS1F_PALINDROME_NULL_DIST.png" not in graphs
    assert "IS1F_CURVATURE.png" not in graphs
    assert "IS1F_MELTING_PROFILE.png" not in graphs


def test_main_analyse_tsds_writes_split_files_and_plots(tmp_path, monkeypatch):
    run_main(monkeypatch, [
        "--mge_name", "IS1F",
        "--summary_file", "data/flanks_ISFinder_IS1F.txt",
        "--flank_len", "50",
        "--outdir", str(tmp_path),
        "--analyse_TSDs",
        "--seed", "1",
    ])

    files = {p.name for p in (tmp_path / "files").glob("*")}
    assert "IS1F_TSD.txt" in files
    assert "IS1F_no_TSD.txt" in files
    assert "IS1F_TSD_SUMMARY.txt" in files

    graphs = {p.name for p in (tmp_path / "graphs").glob("*.png")}
    assert "IS1F_TSD_NULL_DIST.png" in graphs
    assert "IS1F_TSD_LENGTH_DIST.png" in graphs


def test_main_palindrome_centered_uses_flank_boundary_not_hardcoded_50(tmp_path, monkeypatch):
    """--palindrome_centered used to only ever test index 50 as the palindrome
    centre (a bug in find_pals.py), so it silently found nothing unless
    flank_len happened to put the true left/right boundary at 50. The centre
    is now derived from the actual sequence length, so this must find a
    palindrome planted exactly at the boundary of a 20bp flank (boundary at
    index 20 -- nowhere near the old hardcoded 50)."""
    arm = "ACGTAC"
    rc_arm = revcomp(arm)
    left = "G" * 14 + arm       # 20bp, arm ends exactly at the boundary
    right = rc_arm + "G" * 14  # 20bp, arm starts exactly at the boundary

    tsv = tmp_path / "summary.txt"
    write_tsv(tsv, [("acc1", left, right)])

    run_main(monkeypatch, [
        "--mge_name", "TESTPAL",
        "--summary_file", str(tsv),
        "--flank_len", "20",
        "--outdir", str(tmp_path),
        "--analyse_palindromes",
        "--palindrome_centered",
        "--palindrome_min_arm", "3",
        "--palindrome_max_spacer", "0",
    ])

    summary_path = tmp_path / "files" / "TESTPAL_PALINDROME_SUMMARY.txt"
    assert summary_path.exists()
    assert "Observed hits:        1/1 (100.0%)" in summary_path.read_text()


def test_main_seed_is_reproducible_across_runs(tmp_path, monkeypatch):
    """Same --seed should give byte-identical TSD summary files across two
    separate runs (the null-model shuffle is the only source of randomness)."""
    out1, out2 = tmp_path / "run1", tmp_path / "run2"
    common_args = [
        "--mge_name", "IS1F",
        "--summary_file", "data/flanks_ISFinder_IS1F.txt",
        "--flank_len", "50",
        "--analyse_TSDs",
        "--seed", "42",
    ]

    run_main(monkeypatch, common_args + ["--outdir", str(out1)])
    run_main(monkeypatch, common_args + ["--outdir", str(out2)])

    summary1 = (out1 / "files" / "IS1F_TSD_SUMMARY.txt").read_text()
    summary2 = (out2 / "files" / "IS1F_TSD_SUMMARY.txt").read_text()
    assert summary1 == summary2


def test_main_missing_summary_file_exits(tmp_path, monkeypatch):
    with pytest.raises(SystemExit):
        run_main(monkeypatch, [
            "--mge_name", "X",
            "--summary_file", str(tmp_path / "nope.txt"),
            "--flank_len", "10",
            "--outdir", str(tmp_path / "out"),
        ])
