<p align="center">
  <img src="image1350-1.png" width="220" alt="flankalyser logo">
</p>

# flankalyser

Analyse and visualise the DNA sitting either side of an insertion site (e.g. where a transposon/IS element landed). You feed it a table of left/right flanking sequences per accession, it tells you whether anything interesting is going on around the insertion point — repeated sequences (TSDs), inverted repeats (palindromes), unusual base composition, curvature, melting behaviour, that sort of thing.

Still provisional / a bit rough around the edges. Read `python flankalyser.py --help` for the source of truth on flags — this file is the friendlier version.

## Input format

A tab-separated file, one row per accession:

```
acc_name    left_flank    right_flank
```

- `left_flank`/`right_flank` should be the same fixed length (whatever you tell `--flank_len`).
- Multiple flanks for one accession can be pipe-separated (`ACGT...|TTGC...`), and `none` marks a missing side.
- `data/flanks_ISFinder_IS1F.txt` is a real example if you want to see the shape of it.

This is basically the TETyper output format — if your flanks come from somewhere else, just reshape them into this.

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

- **Nucleotide composition** (always runs) — per-position A/C/G/T frequency either side of the insertion site, with a chi-square test against a flat 25% expectation. Asterisks on the plot mark positions where the base composition looks non-random. Worth a look before anything else — if this is flat and boring, the flanks probably aren't telling you much.
- **`--analyse_TSDs`** — looks for target site duplications: a short repeated sequence at the boundary between the left flank and right flank, which is how a lot of transposases mark their landing site. Tests the observed hit rate against a shuffled-pairing null and reports a p-value. A short, consistent TSD length with a low p-value is a decent sign the element actually has one; a flat/noisy distribution suggests either this family doesn't produce TSDs, or the length here isn't detecting them properly.
- **`--analyse_palindromes`** — searches for inverted repeats (a sequence followed by its own reverse complement) near the insertion site, compared to a dinucleotide-shuffled null. Some mobile elements have these as part of their terminal structure. This step is the slow one (repeated shuffling), hence it being optional.
- **`--analyse_bendability`** — estimates DNA backbone curvature/writhe/bend angle at the site using a wedge model. The idea is that some insertion sites might prefer more "bendable" DNA, but this is a physical model, not a measurement, so treat the output as a rough geometric hint rather than a fact.
- **`--analyse_melting`** — looks at how easy the DNA is to pull apart locally: H-bond count and nearest-neighbour stacking energy, position by position, against a shuffled control. A dip in stability right at the insertion boundary could just mean there's more A/T there, or could be reflecting something biologically real — the plot doesn't tell you which, just where.

None of these steps know what they're looking at biologically — they're pattern/statistics tools. Read the numbers, don't over-trust the plots.

## Output

Everything lands under `--outdir`:

```
<outdir>/
  graphs/   *.png, one per analysis, named <MGE_NAME>_<WHAT>.png
  files/    *.txt/.tsv summaries, and the split TSD/palindrome data
```

## Known rough edges

- `--curvature_window` and a few TSD/palindrome parameters are exposed as flags; several others (palindrome arm length, spacer, mismatch tolerance, shuffle counts) are still hardcoded in `flankalyser.py` with a `# TODO` next to them — bump those directly if the defaults don't suit your data.
- GC-content/species-level stuff got pulled out for now; it doesn't need much beyond what's already in the ATB data, just hasn't been rebuilt yet.
- `test/` has skeleton test files (docstrings describing what *should* be tested) rather than full coverage — useful as a checklist, not a safety net yet.
