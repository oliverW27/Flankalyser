<table>
  <tr>
    <td valign="middle"><h1>flankalyser</h1></td>
    <td valign="middle"><img src="image1350-1.png" width="120" alt="flankalyser logo"></td>
  </tr>
</table>

## About

This tool takes in "flanks", i.e. the DNA sequences surrounding some genetic element such as a transposon, gene, or motif. Output includes information and graphs on palindromic content, physical properties (e.g., bendability and curvature), nucleotide and di-nucleotide enrichment, as well as (specific to transposons) target site duplication detection.

## Input format

A tab-separated file, one row per accession:

```
acc_name    left_flank    right_flank
```

- `left_flank`/`right_flank` should be the same fixed length (set by the parameter `--flank_len`).
- `left_flank`/`right_flank` may contain `N` values
- Multiple flanks for one accession can be pipe-separated (`ACGT...|TTGC...`), and `none` marks a missing side.
- See `data_format_example.txt` in `examples/generic`

`helper_scripts` contains a script `convert_to_TETyper_fmt` to convert LexicMap (source here) output into the correct input format

## Requirements

**Flankalyser** requires some basic graphing and calculation libraries, such as numpy, and matplotlib. A conda environment `flankalyser_env` is provided, information on how to use conda is available here (source here).

Run
```bash
conda env create -f environment.yml -n flankalyser_env
conda activate flankalyser_env
```
## Quick start

```bash
python flankalyser.py \
  --mge_name IS1F \
  --summary_file data/flanks_ISFinder_IS1F.txt \
  --flank_len 50 \
  --outdir results/
```

That alone gives you deduplicated flank counts and a nucleotide-composition plot. Everything else is opt-in via flags, because some of it (palindrome shuffling especially) is slow on real datasets:

```bash
python flankalyser.py \
  --mge_name IS1F --summary_file data/flanks_ISFinder_IS1F.txt --flank_len 50 --outdir results/ \
  --analyse_TSDs --analyse_palindromes --analyse_bendability --analyse_melting
```

To run every `flanks_ISFinder_*.txt` file in a directory in one go, there's `run_flankalyser_multi.sh <dir> <flank_len>`.

## What each step does

- **Nucleotide composition** (always runs) — per-position A/C/G/T frequency either side of the insertion site, with a chi-square test against a flat 25% expectation. Asterisks on the plot mark positions where the base composition looks non-random (* : p < 0.05. ** : p < 0.01. *** : p < 0.001). A nucleotide logo is also generated to more easily read any consensus sequence. The number of unique di-nucleotides is compared between the "inner" and "outer" parts of the flanks to see whether there is any enrichment close to the insertion site.

Note - the nucleotide distribution is calculated based on an "expected" baseline of A/C/T/G being equally common. This is not the case in most bacterial genomes. An additional file `gc_content.tsv` can be passed in to use the actual baseline. See (source here) for information on how to generate this file.

- **`--analyse_TSDs`** — looks for target site duplications: a short repeated sequence at the boundary between the left flank and right flank, which is how a lot of transposases mark their landing site. Tests the observed hit rate against a shuffled-pairing null and reports a p-value. A short (ish), consistent TSD length with a low p-value is a decent sign the element actually has one; a flat/noisy distribution suggests either this family doesn't produce TSDs, or they are not presrved in the dataset (e.g., due to post-insertion rearrangements)

- **`--analyse_palindromes`** — searches for palindromes (e.g., AAACGTTT) near the insertion site, compared to a dinucleotide-shuffled null. This parameter may require some fiddling with the parameters `--window` etc.

- **`--analyse_bendability`** — estimates DNA backbone curvature, writhe, and bend angle at the site using a wedge model. The idea is that some transposons might prefer more "bendable" DNA (source here)

- **`--analyse_melting`** — looks at how easy the DNA is to pull apart locally: H-bond count and nearest-neighbour stacking energy, position by position, against a shuffled control.


## Output

Everything lands under `--outdir`:

```
<outdir>/
  graphs/   *.png, one per analysis, named <MGE_NAME>_<WHAT>.png
  files/    *.txt/.tsv summaries, and the split TSD/palindrome data
```

## TODO

- `--curvature_window` and a few TSD/palindrome parameters are exposed as flags; several others (palindrome arm length, spacer, mismatch tolerance, shuffle counts) are still hardcoded in `flankalyser.py` with a `# TODO` next to them — bump those directly if the defaults don't suit your data.
- GC-content/species-level stuff got pulled out for now; it doesn't need much beyond what's already in the ATB data, just hasn't been rebuilt yet.
- `test/` has skeleton test files (docstrings describing what *should* be tested) rather than full coverage — useful as a checklist, not a safety net yet.


## Further analysis, testing, and examples

 - See `helpers` for more scripts that may be helpful, but that require a bit of manual work. For example, species-level analysis as well as more accurate nucleotide distribution using the sequence's actual expected GC content.
- See `tests` for unit tests as well as information on how to run them.
- See `examples` for example input and output data, along with a brief guide on interpretation. 