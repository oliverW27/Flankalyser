"""
Sanity-check script for tsd_analysis.py, run against a real flanks TSV.

Not a behavioral-difference test — this runs the actual pipeline on real data
and prints readable summary statements so you can eyeball them against TSD
lengths you already know to expect (e.g. "IS elements in this family usually
have a 9bp TSD").

Fill in the globals below, then run:
    pytest test_tsd_sanity.py -v -s
"""
import csv
import importlib.util
import sys
import types
from collections import Counter
from pathlib import Path

import pytest

# ---------------------------------------------------------------------------
# FILL THESE IN
# ---------------------------------------------------------------------------
TSD_SCRIPT_PATH = "src/tsd_analysis.py"
FLANK_TSV_PATH = "../flankalyser/IS3/flanks_ISFinder_IS1397.txt"      # acc_name<TAB>left<TAB>right (left/right may be "a|b|c")
FLANK_LEN = 50                             # length of one flank (left or right), before "|" splitting
MGE_NAME = "test_mge"
OUTPUT_DIR = "."
MIN_TSD = 3                                # shortest TSD length to count as a hit
N_SHUFFLES = 100
SEED = 42
EXPECTED_TSD_LEN = None                    # optional: a known/expected TSD length (int) to flag
                                            # against the observed mode, e.g. 9 for many IS elements.
                                            # Leave None to skip that check.
                                            # against the observed mode, e.g. 9 for many IS elements.
                                            # Leave None to skip that check.
# ---------------------------------------------------------------------------


def _load_module(path, name):
    if not Path(path).exists():
        pytest.skip(f"{name} not found at {path!r} — set the path at the top of this file")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def tsd(): # noqa: pytest fixture name matches module purpose, not a redefinition risk here
    # tsd_analysis.py imports `from src.nucl_dist import ...`. If that internal
    # package isn't importable here, stub it so the module still loads — this
    # script only exercises the TSD-detection/significance/plotting path, not
    # the nucleotide-distribution graphs.
    try:
        import src.nucl_dist  # noqa: F401
    except ImportError:
        stub_pkg = types.ModuleType("src")
        stub_mod = types.ModuleType("src.nucl_dist")
        stub_mod.count_nucl = lambda seqs_left, seqs_right, flank_len: ({}, {})
        stub_mod.nucl_dist_stat_sig = lambda counts, gc_content=None: {}
        stub_mod.calc_nucl_props = lambda l, r, flank_len: {n: [0.25] * flank_len for n in "ACGT"}
        stub_mod.adjust_gc_content = lambda props, flank_len: props
        stub_mod.graph = lambda *a, **k: None
        sys.modules["src"] = stub_pkg
        sys.modules["src.nucl_dist"] = stub_mod
    return _load_module(TSD_SCRIPT_PATH, "tsd_analysis")


def _load_flanks(path):
    """Read acc_name/left/right TSV into three parallel lists."""
    acc_names, lefts, rights = [], [], []
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue
            acc_names.append(row[0])
            lefts.append(row[1])
            rights.append(row[2])
    return acc_names, lefts, rights


def _mean(xs):
    return sum(xs) / len(xs) if xs else 0.0


def _median(xs):
    if not xs:
        return 0.0
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


@pytest.mark.skipif(
    not Path(FLANK_TSV_PATH).exists(),
    reason=f"set FLANK_TSV_PATH to a real file (currently {FLANK_TSV_PATH!r})",
)
def test_tsd_sanity_check(tsd):
    acc_names, lefts, rights = _load_flanks(FLANK_TSV_PATH)
    assert acc_names, f"No rows loaded from {FLANK_TSV_PATH}"

    stats = tsd.run_tsd_analysis(
        lefts, rights, acc_names, FLANK_LEN, MGE_NAME, OUTPUT_DIR,
        min_tsd=MIN_TSD, n_shuffles=N_SHUFFLES, seed=SEED,
    )

    total = len(acc_names)
    obs_lengths = stats["obs_lengths"]
    length_counts = Counter(obs_lengths)
    hit_rate = stats["observed_hits"] / total if total else 0.0

    print("\n" + "=" * 60)
    print(f"TSD SANITY CHECK — {MGE_NAME}")
    print("=" * 60)
    print(f"Input file:            {FLANK_TSV_PATH}")
    print(f"Total accessions:      {total}")
    print(f"With a TSD (>= {MIN_TSD}bp): {stats['observed_hits']} ({hit_rate:.1%})")
    print(f"Without a TSD:         {total - stats['observed_hits']} ({1 - hit_rate:.1%})")

    if obs_lengths:
        print(f"\nTSD length — min: {min(obs_lengths)}, max: {max(obs_lengths)}, "
              f"mean: {_mean(obs_lengths):.2f}, median: {_median(obs_lengths):.1f}")
        print("Length breakdown:")
        for length, count in sorted(length_counts.items()):
            pct = 100 * count / len(obs_lengths)
            bar = "#" * max(1, round(pct / 2))
            print(f"  {length:>3}bp: {count:>5} ({pct:5.1f}%)  {bar}")

        most_common_len, most_common_count = length_counts.most_common(1)[0]
        print(f"\nMost common TSD length: {most_common_len}bp "
              f"({most_common_count}/{len(obs_lengths)} hits, "
              f"{100 * most_common_count / len(obs_lengths):.1f}%)")

        if EXPECTED_TSD_LEN is not None:
            if most_common_len == EXPECTED_TSD_LEN:
                print(f"✓ Matches expected TSD length of {EXPECTED_TSD_LEN}bp")
            else:
                print(f"⚠ Expected {EXPECTED_TSD_LEN}bp but most common observed is "
                      f"{most_common_len}bp — worth a closer look")
    else:
        print("\nNo TSDs found at all — check FLANK_LEN, MIN_TSD, and that the")
        print("input file's left/right columns are the sequences you expect.")

    print(f"\nNull mean hits (shuffled pairing): {stats['null_mean']:.2f}/{total}")
    print(f"Enrichment over null:              {stats['enrichment']:.2f}x")
    print(f"p-value (one-sided, {N_SHUFFLES} shuffles):     {stats['p_value']:.4f}")

    out_dir = Path(OUTPUT_DIR)
    print(f"\nFiles written to:  {out_dir.resolve() / 'files'}")
    print(f"Plots written to:  {out_dir.resolve() / 'graphs'}")
    print("=" * 60)

    # Loose sanity assertions — these catch pipeline breakage, not biology.
    assert stats["observed_hits"] >= 0
    assert 0.0 <= stats["p_value"] <= 1.0
