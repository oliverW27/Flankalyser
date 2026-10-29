#!/usr/bin/env python3
"""
Analyses flanking sequences from TETyper output files.
"""

import argparse
import csv
import fcntl
from pathlib import Path
import logging
import sys

# Custom module imports
from src.find_pals import run_palindrome_analysis
from src.tsd_analysis import run_tsd_analysis
from src.nucl_dist import count_nucl, graph, nucl_dist_stat_sig, calc_nucl_props
from src.nucl_enrichment import (
    analyse_dinucleotide_enrichment_weighted)
from src.logo import generate_sequence_logo
from src.dna_bendability import run_bendability_analysis, MODELS as BEND_MODELS
from src.temperature_annealing_calc import run_melting_analysis


def load_gc_content(gc_file, mge_name, logger):
    """Load expected GC content for a specific MGE from a TSV file.

    Args:
        gc_file: Path to TSV with the element name in column 2 and its GC
                 percentage in column 3 (e.g. family, element, mean_gc, ...)
        mge_name: Name of the MGE to look up
        logger: Logger instance

    Returns:
        GC content as decimal (e.g., 0.518 for 51.8%), or None if not found
    """
    if not Path(gc_file).exists():
        logger.error(f"GC file not found: {gc_file}")
        return None

    try:
        with open(gc_file, 'r') as f:
            header = f.readline()  # Skip header
            for line in f:
                cols = line.strip().split('\t')
                if len(cols) >= 3 and cols[1] == mge_name:
                    gc_percent = float(cols[2])
                    gc_decimal = gc_percent / 100.0
                    logger.info(f"Loaded GC content for {mge_name}: {gc_percent}%")
                    return gc_decimal

        logger.warning(f"MGE '{mge_name}' not found in GC file, using default 50%")
        return None
    except Exception as e:
        logger.exception(f"Failed to parse GC file: {e}")
        return None


def log_step(logger, title):
    """Print a banner marking the start of a new pipeline step, to keep console output legible."""
    logger.info("")
    logger.info("=" * 60)
    logger.info(title)
    logger.info("=" * 60)


def parse_summary_file(summary_file, logger):
    """Parse TETyper summary file into separate data lists."""
    acc_names, l_data, r_data = [], [], []

    if not Path(summary_file).exists():
        logger.error(f"Summary file not found: {summary_file}")
        sys.exit(1)

    try:
        with open(summary_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                cols = line.strip().split('\t')
                if len(cols) < 3:
                    logger.warning(
                        f"Skipping malformed line {line_num} in {summary_file}: {line.strip()}"
                    )
                    continue
                acc_names.append(cols[0])
                l_data.append(cols[1])
                r_data.append(cols[2])
    except Exception as e:
        logger.exception(f"Failed while reading summary file: {e}")
        sys.exit(1)

    logger.info(
        f"Parsed {len(acc_names)} entries from {summary_file} "
        f"({len(l_data)} left, {len(r_data)} right flanks)"
    )
    return l_data, r_data, acc_names


def extract_unique_sequences(l_data, r_data, acc_names, flank_len, logger, exclude_n=False):
    """Extract unique sequences from flanking data, excluding 'none' entries."""

    def _contains_n(seq):
        """Check if sequence contains any N characters."""
        return 'N' in seq.upper()

    def _extract_unique(data, names, label):
        unique_seqs, unique_names = [], []
        seq_counts = {}
        skipped_n = 0
        
        for i, sequence in enumerate(data):
            if not sequence:
                logger.debug(f"Skipping empty {label} sequence for {names[i]}")
                continue
            segments = sequence.split('|')
            for segment in segments:
                if segment == 'none':
                    continue
                if len(segment) != flank_len:
                    logger.warning(
                        f"Sequence length mismatch for {names[i]} in {label}: "
                        f"expected {flank_len}, got {len(segment)}"
                    )
                    continue
                if exclude_n and _contains_n(segment):
                    skipped_n += 1
                    continue
                if segment not in unique_seqs:
                    unique_seqs.append(segment)
                    unique_names.append(names[i])
                    seq_counts[segment] = 1
                else:
                    seq_counts[segment] += 1
        
        if exclude_n and skipped_n > 0:
            logger.info(f"Excluded {skipped_n} {label} sequences containing Ns")

        return unique_seqs, unique_names

    l_uniq, l_names = _extract_unique(l_data, acc_names, "left")
    r_uniq, r_names = _extract_unique(r_data, acc_names, "right")

    logger.info(
        f"Extracted {len(l_uniq)} unique left and {len(r_uniq)} unique right sequences"
    )
    return l_uniq, r_uniq, l_names, r_names


def _format_position(v):
    """Render a position with an explicit '+' for positive values (e.g. -7, +8)."""
    return f"+{v}" if v > 0 else str(v)


def _top_dinucleotide(df):
    """Return (dinucleotide, value) for the highest-enrichment cell in a 4x4 ACGT DataFrame."""
    best_val, best_label = None, None
    for a in "ACGT":
        for b in "ACGT":
            val = df.loc[a, b]
            if best_val is None or val > best_val:
                best_val, best_label = val, a + b
    return best_label, best_val


def write_master_table_row(
    master_table_path, mge_name, family, flank_len,
    left_data, right_data, left_seqs, right_seqs,
    left_stats, right_stats, dinuc_dfs, tsd_stats, palindrome_summary,
    logger,
):
    """Append this element's row to a shared master summary table (header
    written once, if the file is new). Every value is derived from data this
    same run already computed in memory -- nothing is re-read from disk.
    """
    # % with a TSD, at whichever length has the highest hit rate
    if tsd_stats and tsd_stats.get('length_counts'):
        total = tsd_stats['total']
        best_len, best_count = max(tsd_stats['length_counts'].items(), key=lambda kv: kv[1])
        tsd_str = f"{best_count / total * 100:.1f}% - {best_len}"
    else:
        tsd_str = "NA"

    # (unique left + unique right) / (total left + total right)
    total_flanks = len(left_data) + len(right_data)
    unique_flanks = len(left_seqs) + len(right_seqs)
    uniqueness_ratio = f"{unique_flanks / total_flanks:.3f}" if total_flanks else "NA"

    # Min/max of all Bonferroni-significant positions, -flank_len..+flank_len
    n_positions = 2 * flank_len
    bonferroni_alpha = 0.05 / n_positions
    sig_positions = [-flank_len + i for i, p in enumerate(left_stats) if p < bonferroni_alpha]
    sig_positions += [i + 1 for i, p in enumerate(right_stats) if p < bonferroni_alpha]
    sig_str = (f"{_format_position(min(sig_positions))}, {_format_position(max(sig_positions))}"
               if sig_positions else "none")

    # % of sequences at/above the palindrome inflection arm length
    pct_infl = palindrome_summary.get('pct_at_inflection') if palindrome_summary else None
    pct_infl_str = f"{pct_infl:.1f}%" if pct_infl is not None else "NA"

    # Top dinucleotide enrichment, left and right
    if dinuc_dfs:
        l_label, l_val = _top_dinucleotide(dinuc_dfs['left'])
        r_label, r_val = _top_dinucleotide(dinuc_dfs['right'])
        dinuc_str = f"{l_label} ({l_val:.2f}) - {r_label} ({r_val:.2f})"
    else:
        dinuc_str = "NA"

    header = ["element", "family", "pct_tsd_best", "uniqueness_ratio",
              "consensus_sequence", "max_consensus_pct", "significant_positions",
              "pct_at_inflection", "top_dinuc_enrichment", "same_pattern_also_in",
              "other_features"]
    row = [mge_name, family, tsd_str, uniqueness_ratio,
           "", "", sig_str, pct_infl_str, dinuc_str, "", ""]

    path = Path(master_table_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Several flankalyser.py runs (different elements/families) may append
    # to this same shared table at once. An exclusive lock around the whole
    # check-then-write sequence below makes that safe: without it, two
    # processes could both see "file doesn't exist yet" and each create it,
    # the second silently wiping out the first's row (or, during a column
    # migration, both read the same N rows and each write back N+1,
    # dropping whichever one loses the race).
    lock_path = path.with_name(path.name + ".lock")
    with open(lock_path, "w") as lockfile:
        fcntl.flock(lockfile, fcntl.LOCK_EX)
        try:
            existing_header = None
            if path.exists():
                with open(path, newline="") as f:
                    existing_header = next(csv.reader(f, delimiter="\t"), [])

            if existing_header is None:
                with open(path, "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(header)
                    writer.writerow(row)
            elif existing_header == header:
                with open(path, "a", newline="") as f:
                    csv.writer(f, delimiter="\t").writerow(row)
            else:
                # The file was created under a different column order (e.g.
                # this script's columns were reordered after that file
                # existed) -- migrate every row by column name rather than
                # desyncing or overwriting it with the new order.
                with open(path, newline="") as f:
                    existing_rows = [dict(r) for r in csv.DictReader(f, delimiter="\t")]
                with open(path, "w", newline="") as f:
                    writer = csv.writer(f, delimiter="\t")
                    writer.writerow(header)
                    for r in existing_rows:
                        writer.writerow([r.get(col, "NA") for col in header])
                    writer.writerow(row)
                logger.warning(f"{path} used a different column order -- migrated "
                                f"{len(existing_rows)} existing rows to the current one")
        finally:
            fcntl.flock(lockfile, fcntl.LOCK_UN)
    logger.info(f"Master table row appended to {path}")


def main():
    """Main pipeline execution."""

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description="Analyse flanking sequences from TETyper output")
    parser.add_argument("--mge_name", required=True,
                        help="Name of the MGE (used in graphs and output files)")
    parser.add_argument("--summary_file",
                        help="Path to TETyper summary file. If not provided, data will be randomly generated.")
    parser.add_argument("--flank_len", required=True, type=int,
                        help="Flanking region length used for TETyper")
    parser.add_argument("--outdir", required=True,
                        help="The directory this program will be run in")
    parser.add_argument("--family", default="",
                        help="ISFinder family this element belongs to (recorded in "
                             "--master_table, if given)")
    parser.add_argument("--master_table",
                        help="Path to a shared TSV summary table -- if given, this run's "
                             "row (TSD rate, uniqueness ratio, significant positions, "
                             "palindrome inflection %%, top dinucleotide enrichment, ...) "
                             "is appended to it (header written once, if the file is new)")
    parser.add_argument("--exclude_n", action="store_true",
                        help="Exclude sequences containing N characters (partial flanks)")
    parser.add_argument("--gc_file",
                        help="TSV file with expected background GC content per element "
                             "(element name in column 2, GC%% in column 3), used to adapt "
                             "the nucleotide-composition significance test and graph "
                             "baseline to the element's actual background instead of a "
                             "flat 25%% per base")
    parser.add_argument("--analyse_TSDs", action="store_true",
                        help="Detect TSDs, split accessions into TSD/no-TSD/trimmed groups, "
                             "write summary files, and generate nucleotide distribution graphs for each group")
    parser.add_argument("--min_tsd", type=int, default=3,
                        help="Minimum TSD length to consider (default: 3)")
    parser.add_argument("--tsd_test_significance", action="store_true",
                        help="Also test the TSD hit rate against a shuffled-pairing null "
                             "and plot it (default: off -- just reports what was found; "
                             "this is the expensive, O(n_seqs)-per-shuffle part)")
    parser.add_argument("--analyse_palindromes", action="store_true",
                        help="Detect inverted repeats/palindromes near insertion sites "
                             "(via src.find_pals) and write a summary + significance plots")
    parser.add_argument("--palindrome_min_arm", type=int, default=6,
                        help="Minimum inverted-repeat arm length (default: 6)")
    parser.add_argument("--palindrome_max_spacer", type=int, default=2,
                        help="Maximum spacer size between inverted-repeat arms (default: 2)")
    parser.add_argument("--palindrome_max_mismatches", type=int, default=0,
                        help="Maximum mismatches allowed per arm (default: 0)")
    parser.add_argument("--palindrome_window_start", type=int, default=None,
                        help="Palindrome search window start, inclusive (default: whole sequence)")
    parser.add_argument("--palindrome_window_end", type=int, default=None,
                        help="Palindrome search window end, exclusive (default: whole sequence)")
    parser.add_argument("--palindrome_n_shuffles", type=int, default=100,
                        help="Number of dinucleotide shuffles for the palindrome null (default: 100)")
    parser.add_argument("--palindrome_min_seq_len", type=int, default=20,
                        help="Skip concatenated left+right sequences shorter than this (default: 20)")
    parser.add_argument("--palindrome_centered", action="store_true",
                        help="Only test the left/right flank boundary as the palindrome centre, "
                             "instead of searching the whole concatenated sequence")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for null-model shuffling (TSD pairing shuffle, "
                             "palindrome dinucleotide shuffle) -- for reproducible runs")
    parser.add_argument("--analyse_bendability", action="store_true",
                        help="Compute the mean local DNA bend profile across insertion "
                             "sites vs. dinucleotide- and mononucleotide-shuffled "
                             "controls (via src.dna_bendability)")
    parser.add_argument("--bend_model", default="nucleosome", choices=BEND_MODELS,
                        help="Bending model used to derive step roll/tilt/twist "
                             "(default: nucleosome)")
    parser.add_argument("--bend_bracket", type=int, default=1,
                        help="Bracket half-width (bp) for the local bend angle: the bend "
                             "at a position is the angle between the base pair normals "
                             "this many steps either side of it (default: 1)")
    parser.add_argument("--bend_n_shuffles", type=int, default=5,
                        help="Shuffles per sequence for each bendability null (default: 5)")
    parser.add_argument("--bend_max_seqs", type=int, default=5000,
                        help="Downsample to at most this many insertion-site sequences "
                             "for the bend profile, 0 for no limit (default: 5000)")
    parser.add_argument("--analyse_melting", action="store_true",
                        help="Compute H-bond count and nearest-neighbor ΔG melting profiles "
                             "vs. a dinucleotide-shuffled control (via src.temperature_annealing_calc)")
    args = parser.parse_args()
    output_dir = args.outdir

    gc_content = None
    if args.gc_file:
        gc_content = load_gc_content(args.gc_file, args.mge_name, logger)

    log_step(logger, "STEP 1: Parsing input")
    logger.info(f"Parsing {args.summary_file}...")
    left_data, right_data, acc_names = parse_summary_file(args.summary_file, logger)

    if not left_data or not right_data:
        logger.error("No valid flank data found. Exiting.")
        sys.exit(1)

    left_seqs, right_seqs, left_names, right_names = extract_unique_sequences(
        left_data, right_data, acc_names, args.flank_len, logger,
        exclude_n=args.exclude_n
    )

    if not left_seqs and not right_seqs:
        logger.error("No unique flank sequences found. Exiting.")
        sys.exit(1)

    MIN_UNIQUE_FLANKS = 1
    if len(left_seqs) < MIN_UNIQUE_FLANKS or len(right_seqs) < MIN_UNIQUE_FLANKS:
        msg = (f"SKIPPED: {args.mge_name} has only {len(left_seqs)} unique left / "
               f"{len(right_seqs)} unique right flanks (minimum {MIN_UNIQUE_FLANKS} "
               f"required) -- no results generated.")
        logger.warning(msg)
        print(msg)
        sys.exit(0)

    try:
        Path(f'{output_dir}/graphs').mkdir(parents=True, exist_ok=True)
        Path(f'{output_dir}/files').mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created output directories under {output_dir}")
    except Exception as e:
        logger.exception(f"Could not create output directories: {e}")
        sys.exit(1)

    counts_path = Path(output_dir) / "files" / f"{args.mge_name}_FLANK_COUNTS.tsv"
    with open(counts_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"left_total\t{len(left_data)}\n")
        f.write(f"left_unique\t{len(left_seqs)}\n")
        f.write(f"right_total\t{len(right_data)}\n")
        f.write(f"right_unique\t{len(right_seqs)}\n")
    logger.info(f"Flank counts written to {counts_path}")

    log_step(logger, "STEP 2: Nucleotide composition analysis")
    try:
        left_counts, right_counts = count_nucl(left_seqs, right_seqs, args.flank_len)
        left_stats = nucl_dist_stat_sig(left_counts, gc_content=gc_content)
        right_stats = nucl_dist_stat_sig(right_counts, gc_content=gc_content)
        logger.info("Calculated nucleotide distributions successfully.")
        combined_props = calc_nucl_props(left_counts, right_counts, args.flank_len)

        # Per-position significance (Bonferroni-corrected chi-square), written
        # out so downstream reporting can find which positions (relative to
        # the insertion site, -flank_len..+flank_len) are significant without
        # redoing the stats test.
        n_positions = 2 * args.flank_len
        bonferroni_alpha = 0.05 / n_positions
        sig_path = Path(output_dir) / "files" / f"{args.mge_name}_SIGNIFICANT_POSITIONS.tsv"
        with open(sig_path, "w") as f:
            f.write("position\tp_value\tsignificant\n")
            for i, p in enumerate(left_stats):
                position = -args.flank_len + i
                f.write(f"{position}\t{p}\t{p < bonferroni_alpha}\n")
            for i, p in enumerate(right_stats):
                position = i + 1
                f.write(f"{position}\t{p}\t{p < bonferroni_alpha}\n")
        logger.info(f"Significant positions written to {sig_path}")

    except Exception as e:
        logger.exception(f"Error calculating nucleotide distributions: {e}")
        sys.exit(1)

    graph(
        args.flank_len,
        args.mge_name,
        combined_props,
        left_stats,
        right_stats,
        outdir=args.outdir,
        total_left=len(left_data),
        unique_left=len(left_seqs),
        total_right=len(right_data),
        unique_right=len(right_seqs),
        gc_content=gc_content,
    )

    log_step(logger, "STEP 3: Dinucleotide enrichment & sequence logo")
    dinuc_dfs = analyse_dinucleotide_enrichment_weighted(left_seqs, right_seqs, output_dir, args.flank_len, args.mge_name)
    try:
        generate_sequence_logo(left_seqs, right_seqs, output_dir, args.mge_name)
    except Exception as e:
        logger.exception(f"Sequence logo generation failed: {e}")

    tsd_stats = None
    if args.analyse_TSDs:
        log_step(logger, "STEP 4: TSD analysis")
        try:
            tsd_stats = run_tsd_analysis(
                left_data,
                right_data,
                acc_names,
                args.flank_len,
                args.mge_name,
                output_dir,
                min_tsd=args.min_tsd,
                seed=args.seed,
                test_significance=args.tsd_test_significance,
            )
        except Exception as e:
            logger.exception(f"TSD analysis failed: {e}")

    palindrome_summary = None
    if args.analyse_palindromes:
        log_step(logger, "STEP 5: Palindrome / inverted-repeat analysis")
        try:
          palindrome_summary = run_palindrome_analysis(
            input_tsv=args.summary_file,
            mge_name=args.mge_name,
            output_dir=output_dir,
            min_arm=args.palindrome_min_arm,
            max_arm=25,  # new — add to argparse if you want it configurable
            max_spacer=args.palindrome_max_spacer,
            max_mismatches=args.palindrome_max_mismatches,
            n_shuffles=args.palindrome_n_shuffles,
            threshold=2.0,  # new — add to argparse if you want it configurable
            seed=args.seed,
            min_seq_len=args.palindrome_min_seq_len,
            centered=args.palindrome_centered,
        )
        except Exception as e:
            logger.exception(f"Palindrome analysis failed: {e}")

    if args.analyse_bendability:
        log_step(logger, "STEP 6: DNA bendability analysis")
        try:
            run_bendability_analysis(
                left_data, right_data, acc_names,
                args.flank_len, args.mge_name, output_dir,
                model=args.bend_model,
                n_b=args.bend_bracket,
                n_shuffles=args.bend_n_shuffles,
                max_seqs=args.bend_max_seqs,
                seed=args.seed,
            )
        except Exception as e:
            logger.exception(f"Bendability analysis failed: {e}")

    if args.analyse_melting:
        log_step(logger, "STEP 7: Melting / annealing analysis")
        try:
            run_melting_analysis(left_seqs, right_seqs, args.mge_name, output_dir)
        except Exception as e:
            logger.exception(f"Melting analysis failed: {e}")

    if args.master_table:
        log_step(logger, "STEP 8: Master table row")
        try:
            write_master_table_row(
                args.master_table, args.mge_name, args.family, args.flank_len,
                left_data, right_data, left_seqs, right_seqs,
                left_stats, right_stats, dinuc_dfs, tsd_stats, palindrome_summary,
                logger,
            )
        except Exception as e:
            logger.exception(f"Failed to write master table row: {e}")

    log_step(logger, "DONE")


if __name__ == "__main__":
    main()