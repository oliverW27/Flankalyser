#!/usr/bin/env python3
"""
Analyses flanking sequences from TETyper output files.
"""

import argparse
import random
from pathlib import Path
import logging
import sys

# Custom module imports
from src.kmer_count import plot_unique_kmer_counts_both
from src.palindrome_analysis import process_summary_file_centered_palindromes
from src.tsd_detection import detect_tsd_all_vs_all
from src.tsd_analysis import run_tsd_analysis
from src.nucl_dist import count_nucl, graph, nucl_dist_stat_sig, calc_nucl_props, adjust_gc_content
from src.nucl_enrichment import (
    analyse_dinucleotide_enrichment_weighted,
    analyse_trinucleotide_enrichment_weighted
)
from src.logo import generate_sequence_logo


def load_gc_content(gc_file, mge_name, logger):
    """Load GC content for a specific MGE from TSV file.
    
    Args:
        gc_file: Path to TSV with columns 'element' and 'expected_gc'
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
                if len(cols) >= 2 and cols[1] == mge_name:
                    gc_percent = float(cols[2])
                    gc_decimal = gc_percent / 100.0
                    logger.info(f"Loaded GC content for {mge_name}: {gc_percent}%")
                    return gc_decimal
        
        logger.warning(f"MGE '{mge_name}' not found in GC file, using default 50%")
        return None
    except Exception as e:
        logger.exception(f"Failed to parse GC file: {e}")
        return None


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


def extract_unique_sequences(l_data, r_data, acc_names, flank_len, logger, 
                              remove_likely_snps=False, exclude_n=False):
    """Extract unique sequences from flanking data, excluding 'none' entries."""

    def _hamming_distance(seq1, seq2):
        """Calculate Hamming distance between two sequences of equal length."""
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

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
        
        if remove_likely_snps:
            unique_seqs, unique_names = _remove_snp_duplicates(
                unique_seqs, unique_names, seq_counts, label, logger
            )
        
        return unique_seqs, unique_names

    def _remove_snp_duplicates(seqs, names, counts, label, logger):
        """Remove sequences that differ by exactly one SNP, keeping the most frequent."""
        merged = {}
        final_seqs = []
        final_names = []
        
        for i, seq in enumerate(seqs):
            if i in merged:
                continue
            
            similar_indices = [i]
            for j in range(i + 1, len(seqs)):
                if j not in merged and _hamming_distance(seqs[i], seqs[j]) == 1:
                    similar_indices.append(j)
            
            best_idx = max(similar_indices, key=lambda idx: counts[seqs[idx]])
            
            if len(similar_indices) > 1:
                similar_seqs = [seqs[idx] for idx in similar_indices]
                logger.info(
                    f"Merging {len(similar_seqs)} similar {label} sequences "
                    f"(differ by 1 SNP): {similar_seqs} -> keeping '{seqs[best_idx]}' "
                    f"(count: {counts[seqs[best_idx]]})"
                )
            
            for idx in similar_indices:
                merged[idx] = best_idx
            
            final_seqs.append(seqs[best_idx])
            final_names.append(names[best_idx])
        
        return final_seqs, final_names

    l_uniq, l_names = _extract_unique(l_data, acc_names, "left")
    r_uniq, r_names = _extract_unique(r_data, acc_names, "right")

    logger.info(
        f"Extracted {len(l_uniq)} unique left and {len(r_uniq)} unique right sequences"
    )
    return l_uniq, r_uniq, l_names, r_names


def generate_test_data(flank_len, num_accs, outfile, logger):
    """Generate random test data for pipeline testing."""
    nucleotides = ['A', 'G', 'C', 'T']

    def _make_random_flank():
        seq = ''.join(random.choices(nucleotides, k=flank_len))
        return seq

    try:
        accessions = []
        for i in range(num_accs):
            acc_name = f"acc_{i}"
            left_flank = _make_random_flank()
            right_flank = _make_random_flank()
            accessions.append(f"{acc_name}\t{left_flank}\t{right_flank}")

        with open(outfile, "w") as f:
            for acc in accessions:
                f.write(f"{acc}\n")

        logger.info(f"Generated {num_accs} test accessions in {outfile}")
    except Exception as e:
        logger.exception(f"Failed to generate test data: {e}")
        sys.exit(1)


def run_analysis_suite(summary_file, left_seqs, right_seqs, output_dir, flank_len, logger, combined_props=None):
    """Run comprehensive analysis suite on flanking sequences."""
    if combined_props is None:
        combined_props = {}

    analyses = [
        lambda: analyse_dinucleotide_enrichment_weighted(left_seqs, right_seqs, combined_props, output_dir, flank_len),
        lambda: analyse_trinucleotide_enrichment_weighted(left_seqs, right_seqs, combined_props, output_dir, flank_len),
        lambda: generate_sequence_logo(left_seqs, right_seqs, output_dir),
        lambda: process_summary_file_centered_palindromes(left_seqs, right_seqs, combined_props, summary_file),
        lambda: plot_unique_kmer_counts_both(left_seqs, right_seqs, combined_props, flank_len),
        lambda: detect_tsd_all_vs_all(flank_len, combined_props, left_seqs, right_seqs, output_dir),
    ]

    for i, analysis in enumerate(analyses, 1):
        try:
            logger.debug(f"Running analysis step {i}")
            analysis()
        except Exception as e:
            logger.exception(f"Analysis step {i} failed: {e}")


def create_default_proportions(flank_len):
    """Create default nucleotide proportions (equal distribution)."""
    default_prop = [0.25] * flank_len + ["nan"] + [0.25] * flank_len
    return {"A": default_prop, "C": default_prop, "G": default_prop, "T": default_prop}


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
    parser.add_argument("--test", action="store_true",
                        help="Run pipeline on randomly generated test data")
    parser.add_argument("--adjust_proportions", action="store_true",
                        help="Use calculated nucleotide proportions (default: equal)")
    parser.add_argument("--adjust_GC", action="store_true",
                        help="Adjust for background GC content using outer 10bp baseline")
    parser.add_argument("--outdir", required=True,
                        help="The directory this program will be run in")
    parser.add_argument("--rm_likely_snps", action="store_true",
                        help="De-duplicate flanks with only one SNP")
    parser.add_argument("--exclude_n", action="store_true",
                        help="Exclude sequences containing N characters (partial flanks)")
    parser.add_argument("--gc_file",
                        help="TSV file with expected GC content per element (columns: element, expected_gc)")
    parser.add_argument("--analyse_TSDs", action="store_true",
                        help="Detect TSDs, split accessions into TSD/no-TSD/trimmed groups, "
                             "write summary files, and generate nucleotide distribution graphs for each group")
    parser.add_argument("--min_tsd", type=int, default=3,
                        help="Minimum TSD length to consider (default: 3)")
    args = parser.parse_args()
    output_dir = args.outdir

    try:
        Path(f'{output_dir}/graphs').mkdir(parents=True, exist_ok=True)
        Path(f'{output_dir}/files').mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created output directories under {output_dir}")
    except Exception as e:
        logger.exception(f"Could not create output directories: {e}")
        sys.exit(1)

    # Load GC content if provided
    gc_content = None
    if args.gc_file:
        gc_content = load_gc_content(args.gc_file, args.mge_name, logger)

    # Handle test mode
    if args.test or not args.summary_file:
        logger.info("Generating test data...")
        generate_test_data(args.flank_len, 1000, "test_summary.txt", logger)
        args.summary_file = "test_summary.txt"

    # Parse input data
    logger.info(f"Parsing {args.summary_file}...")
    left_data, right_data, acc_names = parse_summary_file(args.summary_file, logger)

    if not left_data or not right_data:
        logger.error("No valid flank data found. Exiting.")
        sys.exit(1)

    left_seqs, right_seqs, left_names, right_names = extract_unique_sequences(
        left_data, right_data, acc_names, args.flank_len, logger,
        remove_likely_snps=args.rm_likely_snps,
        exclude_n=args.exclude_n
    )

    if not left_seqs and not right_seqs:
        logger.error("No unique flank sequences found. Exiting.")
        sys.exit(1)

    # Analyse nucleotide composition
    try:
        left_counts, right_counts = count_nucl(left_seqs, right_seqs, args.flank_len)
        left_stats = nucl_dist_stat_sig(left_counts, gc_content=gc_content)
        right_stats = nucl_dist_stat_sig(right_counts, gc_content=gc_content)
        logger.info("Calculated nucleotide distributions successfully.")
    except Exception as e:
        logger.exception(f"Error calculating nucleotide distributions: {e}")
        sys.exit(1)

    # Calculate or use default proportions
    try:
        if args.adjust_proportions:
            combined_props = calc_nucl_props(left_counts, right_counts, args.flank_len)
            logger.info("Using adjusted nucleotide proportions.")
        else:
            combined_props = create_default_proportions(args.flank_len)
            logger.info("Using default equal nucleotide proportions.")
    except Exception as e:
        logger.exception(f"Error creating nucleotide proportions: {e}")
        sys.exit(1)

    # Apply GC adjustment if requested
    if args.adjust_GC:
        try:
            combined_props = adjust_gc_content(combined_props, args.flank_len)
            logger.info("Applied GC content adjustment.")
        except Exception as e:
            logger.exception(f"Error applying GC adjustment: {e}")
            sys.exit(1)

    # Generate nucleotide distribution graph
    try:
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
        logger.info("Nucleotide distribution graph generated.")
    except Exception as e:
        logger.exception(f"Failed to generate graph: {e}")

    # Run extended analysis
    #run_analysis_suite(args.summary_file, left_seqs, right_seqs, output_dir, args.flank_len, logger, combined_props)

    # TSD analysis
    if args.analyse_TSDs:
        logger.info("Running TSD analysis...")
        try:
            run_tsd_analysis(
                left_data,
                right_data,
                acc_names,
                args.flank_len,
                args.mge_name,
                output_dir,
                logger,
                adjust_proportions=args.adjust_proportions,
                adjust_gc=args.adjust_GC,
                gc_content=gc_content,
                min_tsd=args.min_tsd,
            )
        except Exception as e:
            logger.exception(f"TSD analysis failed: {e}")

    logger.info("Analysis complete!")


if __name__ == "__main__":
    main()