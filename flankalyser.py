#!/usr/bin/env python3
"""
Analyses flanking sequences from TETyper output files.
"""

import argparse
from datetime import datetime
from importlib import metadata
from pathlib import Path
import logging
import shlex
import sys
import matplotlib
matplotlib.use("Agg")

# Custom module imports
from src.find_pals import run_palindrome_analysis, MAX_ARM
from src.tsd_analysis import run_tsd_analysis
from src.nucl_dist import (count_nucl, graph_nucl_dist, nucl_dist_stat_sig, calc_nucl_props,
                           expected_proportions)
from src.nucl_enrichment import (
    analyse_dinucleotide_enrichment_weighted)
from src.logo import generate_sequence_logo
from src.dna_bendability import run_bendability_analysis
from src.temperature_annealing_calc import run_melting_analysis
from src.utils import parse_summary_file

MIN_UNIQUE_FLANKS = 100


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
            next(f, None)  # Skip header
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
    """Print a banner marking the start of a new pipeline step, to keep console output somewhat legible"""
    logger.info("")
    logger.info("-" * 50)
    logger.info(title)
    logger.info("-" * 50)


def extract_unique_sequences(l_data, r_data, acc_names, flank_len, logger, exclude_n=False):
    """Extract unique flank sequences of the expected length (optionally excluding any with N)."""

    def _contains_n(seq):
        """Check if sequence contains any N characters."""
        return 'N' in seq.upper()

    def _extract_unique(data, names, label):
        unique_seqs = []
        seen = set()
        skipped_n = 0
        
        for i, sequence in enumerate(data):
            if not sequence:
                logger.debug(f"Skipping empty {label} sequence for {names[i]}")
                continue
            segments = sequence.split('|')
            for segment in segments:
                if len(segment) != flank_len:
                    logger.warning(
                        f"Sequence length mismatch for {names[i]} in {label}: "
                        f"expected {flank_len}, got {len(segment)}"
                    )
                    continue
                if exclude_n and _contains_n(segment):
                    skipped_n += 1
                    continue
                if segment not in seen:
                    seen.add(segment)
                    unique_seqs.append(segment)
        
        if exclude_n and skipped_n > 0:
            logger.info(f"Excluded {skipped_n} {label} sequences containing Ns")

        return unique_seqs

    l_uniq = _extract_unique(l_data, acc_names, "left")
    r_uniq = _extract_unique(r_data, acc_names, "right")

    logger.info(
        f"Extracted {len(l_uniq)} unique left and {len(r_uniq)} unique right sequences"
    )
    return l_uniq, r_uniq


def setup_logging(outdir, mge_name):
    """Send log messages to the screen and to <outdir>/<mge_name>_<YYYYmmdd-HHMMSS>.log.

    The file copy adds a timestamp to each line; the screen output is unchanged.
    Returns the logger and the path to the log file.
    """
    Path(outdir).mkdir(parents=True, exist_ok=True)
    log_path = Path(outdir) / f"{mge_name}_{datetime.now():%Y%m%d-%H%M%S}.log"

    screen = logging.StreamHandler()
    screen.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    file = logging.FileHandler(log_path)
    file.setFormatter(logging.Formatter('%(asctime)s %(levelname)s: %(message)s',
                                        datefmt='%Y-%m-%d %H:%M:%S'))
    logging.basicConfig(level=logging.INFO, handlers=[screen, file])
    logger = logging.getLogger(__name__)

    # Python warnings (e.g. from numpy or matplotlib) normally print straight to the screen;
    # send them through logging instead so they also reach the log file
    logging.captureWarnings(True)

    # Same for crashes that nothing else catches: log the traceback rather than only printing it
    def log_uncaught(exc_type, exc_value, exc_tb):
        if issubclass(exc_type, KeyboardInterrupt):
            logger.error("Interrupted by user")
            return
        logger.critical("Unexpected error -- run stopped",
                        exc_info=(exc_type, exc_value, exc_tb))
    sys.excepthook = log_uncaught

    return logger, log_path


def log_run_info(logger, log_path):
    """Record the command, Python and package versions at the top of the log."""
    versions = []
    for pkg in ("numpy", "pandas", "scipy", "matplotlib", "seaborn", "logomaker"):
        try:
            versions.append(f"{pkg} {metadata.version(pkg)}")
        except metadata.PackageNotFoundError:
            versions.append(f"{pkg} not found")
    logger.info(f"Command: {shlex.join(sys.argv)}")
    logger.info(f"Python {sys.version.split()[0]}; " + ", ".join(versions))
    logger.info(f"Log file: {log_path}")


def main():
    """Main pipeline execution."""

    parser = argparse.ArgumentParser(description="Analyse flanking sequences of MGEs")
    parser.add_argument("--mge_name", required=True,
                        help="Name of the element; used in graph titles, output file names "
                             "and the log file name, and to look up its GC content in --gc_file")
    parser.add_argument("--summary_file", required=True,
                        help="Tab-separated input file, one line per genome: name, left flank, "
                             "right flank. A field may hold several flanks separated by '|'. "
                             "Missing sequence should be filled with N to the full flank length. "
                             "Bases are converted to uppercase, and any character other than "
                             "A/C/G/T/N is treated as N")
    parser.add_argument("--flank_len", required=True, type=int,
                        help="Length of each flank in bp; flanks of any other length are skipped")
    parser.add_argument("--outdir", required=True,
                        help="Output directory (created if needed). Results go in graphs/ and "
                             "files/, and the log file goes in the directory itself")
    parser.add_argument("--exclude_n", action="store_true",
                        help="Leave out flanks containing N in the nucleotide-composition "
                             "and sequence-logo steps, and in the 100-unique-flank minimum "
                             "(the other analyses always leave out flanks, or TSDs, containing N)")
    parser.add_argument("--gc_file",
                        help="TSV file with a header line and the expected background GC content "
                             "per element (element name in column 2, GC%% in column 3). Sets the "
                             "expected base proportions for the nucleotide-composition test and "
                             "graph, instead of a flat 25%% per base. Elements not in the file "
                             "fall back to 25%% per base")
    parser.add_argument("--analyse_TSDs", action="store_true",
                        help="Find the longest target site duplication at each unique left/right "
                             "flank pair and write a TSD summary")
    parser.add_argument("--min_tsd", type=int, default=3,
                        help="Minimum TSD length to count, with --analyse_TSDs "
                             "(default: 3; maximum is fixed at 25)")
    parser.add_argument("--tsd_test_significance", action="store_true",
                        help="With --analyse_TSDs, also test the TSD hit rate against the same "
                             "flanks with the boundary shifted 10-25 bp, add the result to the "
                             "TSD summary, and plot TSD lengths observed vs. shifted")
    parser.add_argument("--analyse_palindromes", action="store_true",
                        help="Join each unique left/right flank pair, find the longest inverted "
                             "repeat in it, and compare against dinucleotide-shuffled sequences; "
                             "writes the enrichment curve and a plot. Pairs containing N are skipped")
    parser.add_argument("--palindrome_min_arm", type=int, default=3,
                        help="Shortest inverted-repeat arm in bp to count as a palindrome, and the "
                             "first arm length on the enrichment curve, with --analyse_palindromes "
                             "(default: 3)")
    parser.add_argument("--palindrome_max_spacer", type=int, default=2,
                        help="Maximum spacer between the two inverted-repeat arms in bp, with "
                             "--analyse_palindromes (default: 2)")
    parser.add_argument("--palindrome_max_mismatches", type=int, default=0,
                        help="Maximum non-complementary base pairs allowed between the two arms "
                             "of an inverted repeat, with --analyse_palindromes (default: 0)")
    parser.add_argument("--analyse_bendability", action="store_true",
                        help="Plot the mean local DNA bend profile (Goodsell & Dickerson "
                             "nucleosome model) of each flank vs. a dinucleotide-shuffled control. "
                             "Flanks containing N are skipped; at most 5000 randomly chosen "
                             "flanks per side are used")
    parser.add_argument("--analyse_melting", action="store_true",
                        help="Plot the nearest-neighbour stacking ΔG profile (SantaLucia 1998) "
                             "of each flank vs. a dinucleotide-shuffled control. "
                             "Flanks containing N are skipped")
    args = parser.parse_args()
    if not 1 <= args.palindrome_min_arm <= MAX_ARM:
        parser.error(f"--palindrome_min_arm must be between 1 and {MAX_ARM}")
    output_dir = args.outdir

    logger, log_path = setup_logging(output_dir, args.mge_name)
    log_run_info(logger, log_path)

    # Read in files
    log_step(logger, "STEP 1: Parsing input")
    gc_content = None
    if args.gc_file:
        gc_content = load_gc_content(args.gc_file, args.mge_name, logger)

    logger.info(f"Parsing {args.summary_file}...")
    left_data, right_data, acc_names = parse_summary_file(args.summary_file)

    if not left_data or not right_data:
        logger.error("No valid flank data found. Exiting.")
        sys.exit(1)

    left_seqs, right_seqs = extract_unique_sequences(
        left_data, right_data, acc_names, args.flank_len, logger,
        exclude_n=args.exclude_n
    )

    if not left_seqs and not right_seqs:
        logger.error("No unique flank sequences found. Exiting.")
        sys.exit(1)

    if len(left_seqs) < MIN_UNIQUE_FLANKS or len(right_seqs) < MIN_UNIQUE_FLANKS:
        msg = (f"SKIPPED: {args.mge_name} has only {len(left_seqs)} unique left / "
               f"{len(right_seqs)} unique right flanks (minimum {MIN_UNIQUE_FLANKS} "
               f"required) -- no results generated.")
        logger.warning(msg)
        sys.exit(0)

    try:
        Path(f'{output_dir}/graphs').mkdir(parents=True, exist_ok=True)
        Path(f'{output_dir}/files').mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created output directories under {output_dir}")
    except Exception as e:
        logger.exception(f"Could not create output directories: {e}")
        sys.exit(1)

    # total: every left/right flank pair in the input, with pipe-separated pairs counted separately
    # unique_pairs: distinct left+right pairs of the right length (what the TSD step uses)
    total_pairs = sum(min(len(l.split('|')), len(r.split('|'))) for l, r in zip(left_data, right_data))
    n_unique_pairs = len({(l, r) for lf, rf in zip(left_data, right_data)
                          for l, r in zip(lf.upper().split('|'), rf.upper().split('|'))
                          if len(l) == args.flank_len and len(r) == args.flank_len})
    counts_path = Path(output_dir) / "files" / f"{args.mge_name}_FLANK_COUNTS.tsv"
    with open(counts_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total\t{total_pairs}\n")
        f.write(f"left_unique\t{len(left_seqs)}\n")
        f.write(f"right_unique\t{len(right_seqs)}\n")
        f.write(f"unique_pairs\t{n_unique_pairs}\n")
    logger.info(f"Flank counts written to {counts_path}")

    log_step(logger, "STEP 2: Nucleotide composition analysis")
    try:
        left_counts, right_counts = count_nucl(left_seqs, right_seqs, args.flank_len)
        left_stats = nucl_dist_stat_sig(left_counts, gc_content=gc_content)
        right_stats = nucl_dist_stat_sig(right_counts, gc_content=gc_content)
        logger.info("Calculated nucleotide distributions successfully.")
        combined_props = calc_nucl_props(left_counts, right_counts)
        n_positions = 2 * args.flank_len
        bonferroni_alpha = 0.05 / n_positions
        expected = expected_proportions(gc_content)
        if gc_content is not None:
            expected_note = (f"from GC file {args.gc_file} ({gc_content:.1%} GC: "
                             f"A/T {expected['A']:.1%}, C/G {expected['C']:.1%} each)")
        elif args.gc_file:
            expected_note = f"flat 25% per base ({args.mge_name} not found in GC file {args.gc_file})"
        else:
            expected_note = "flat 25% per base (no GC file given)"

        sig_path = Path(output_dir) / "files" / f"{args.mge_name}_SIGNIFICANT_POSITIONS.tsv"
        with open(sig_path, "w") as f:
            f.write(f"# Expected base proportions: {expected_note}\n")
            f.write("position\tp_value\tsignificant\t"
                    + "\t".join(f"{b}_minus_expected" for b in "ACGT") + "\n")
            for counts, stats, positions in ((left_counts, left_stats, range(-args.flank_len, 0)),
                                             (right_counts, right_stats, range(1, args.flank_len + 1))):
                for i, (position, p) in enumerate(zip(positions, stats)):
                    # Observed minus expected proportion of each base (+ = over-represented)
                    n = sum(counts[b][i] for b in "ACGT")
                    diffs = [f"{counts[b][i] / n - expected[b]:+.4f}" if n else "NA" for b in "ACGT"]
                    f.write(f"{position}\t{p}\t{p < bonferroni_alpha}\t" + "\t".join(diffs) + "\n")
        logger.info(f"Significant positions written to {sig_path}")

    except Exception as e:
        logger.exception(f"Error calculating nucleotide distributions: {e}")
        sys.exit(1)

    graph_nucl_dist(
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
    analyse_dinucleotide_enrichment_weighted(left_seqs, right_seqs, output_dir, args.mge_name)
    try:
        generate_sequence_logo(left_counts, right_counts, output_dir, args.mge_name)
    except Exception as e:
        logger.exception(f"Sequence logo generation failed: {e}")

    if args.analyse_TSDs:
        log_step(logger, "STEP 4: TSD analysis")
        try:
            run_tsd_analysis(
                left_data,
                right_data,
                acc_names,
                args.flank_len,
                args.mge_name,
                output_dir,
                min_tsd=args.min_tsd,
                test_significance=args.tsd_test_significance,
            )
        except Exception as e:
            logger.exception(f"TSD analysis failed: {e}")

    if args.analyse_palindromes:
        log_step(logger, "STEP 5: Palindrome / inverted-repeat analysis")
        try:
          run_palindrome_analysis(
            left_data, right_data, args.flank_len,
            mge_name=args.mge_name,
            output_dir=output_dir,
            min_arm=args.palindrome_min_arm,
            max_spacer=args.palindrome_max_spacer,
            max_mismatches=args.palindrome_max_mismatches,
            threshold=2.0,
        )
        except Exception as e:
            logger.exception(f"Palindrome analysis failed: {e}")

    if args.analyse_bendability:
        log_step(logger, "STEP 6: DNA bendability analysis")
        try:
            run_bendability_analysis(
                left_seqs, right_seqs, args.flank_len, args.mge_name, output_dir,
            )
        except Exception as e:
            logger.exception(f"Bendability analysis failed: {e}")

    if args.analyse_melting:
        log_step(logger, "STEP 7: Melting / annealing analysis")
        try:
            run_melting_analysis(left_seqs, right_seqs, args.mge_name, output_dir)
        except Exception as e:
            logger.exception(f"Melting analysis failed: {e}")

    log_step(logger, "DONE")


if __name__ == "__main__":
    main()