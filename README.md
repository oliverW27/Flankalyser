<p align="center">
  <img src="logo.png" width="180" alt="flankalyser logo">
</p>
<h1 align="center">Flankalyser</h1>
<p align="center">Analysis of the DNA flanking bacterial transposable elements</p>

## About

This tool takes in "flanks", i.e. the DNA sequences surrounding some genetic element such as a transposon, gene, or motif. Output includes information and graphs on palindromic content, physical properties (e.g., bendability), nucleotide and dinucleotide enrichment, as well as (specific to transposons) target site duplication detection.

## Dependencies
This tool was tested with Python 3.12 and requires:
  - numpy
  - pandas
  - scipy
  - matplotlib
  - seaborn
  - logomaker

A [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) .yml file is provided. To use it:
```bash
conda env create -f env.yml
conda activate flankalyser
```
## Input
The input file must contain the flanking contexts of a given MGE in the following, tab-delimated format:

```tsv
[unique ID]   [left flank 1]|[optional: left flank 2]|...   [right flank 1]|[optional: right flank 2]...
```
All flanks must be of the same length, but may be padded with `N`. Note the `--flank-len` argument must be accurate for the actual length of the flanks, else all input will be skipped.

**Valid input** (N and invalid characters are ignored):
```tsv
accession_1    NNAGCT    ACCNNN
accession_2    NNNNNN    AGTTTT
accession_3    AGNCCC    TTAANA
accession_4    HELLO1    world2
accession_5    agcttt    tgcnnn
```
- `accession_4`: characters not in [A,G,C,T] are converted to N
- `accession_5`: lower-case nucleotides are fine!

**Problem input** (logs a warning):
```tsv
accession_1    NNAGCT
accession_2    ACG    AGTTTT
accession_3    AGNCCC    AGCTTA|ATTGGG
```
- `accession_1`: only one flank column, so the whole row is skipped.
- `accession_2`: the left flank is the wrong length, so it is skipped. The right flank is still used in the analyses that look at each flank on its own (nucleotide composition, sequence logo, dinucleotide enrichment, bendability, melting), but not in the TSD or palindrome analyses
- `accession_3`: different numbers of left and right flanks. All three flanks are used in the analyses that look at each flank on its own. The TSD and palindrome analyses pair the flanks in order (`AGNCCC` with `AGCTTA`) and leave out the unpaired `ATTGGG`.

The expected way of data generation for this tool is using [LexicMap](https://github.com/shenwei356/LexicMap) to query a database of assemblies, such as [AllTheBacteria](https://allthebacteria.org/).
A script `convert_lexicmap_fmt` is provided in `helpers/` to convert the output format to one compatible with Flankalyser.

### Example input generation workflow
```bash
# Query
lexicmap search \
	-d atb_dir/atb.lmi \
	MGE.fasta \
	-o output.fasta.lexicmap.tsv \
	--align-min-match-pident 80 \
	--min-qcov-per-hsp 70 \
	--min-qcov-per-genome 70

# Extract flanking contexts
lexicmap utils subseq -e -d atb_dir/atb.lmi \
  -f output.fasta.lexicmap.tsv \
  -o MGE_flanks.fasta \
  -U 50 -D 50 # For flank length 50

# Convert to Flankalyser expected format
python3 convert_lexicmap_fmt.py --tsv output.fasta.lexicmap.tsv \
 --input MGE_flanks.fasta \
 --output_dir flanks \
 --flank_len 50

```
## Options

| Option | Type | Required | Default |
|---|---|---|---|
| `--mge_name` | string | yes | – |
| `--summary_file` | path | yes | – |
| `--flank_len` | integer | yes | – |
| `--outdir` | path | yes | – |
| `--exclude_n` | flag | no | off |
| `--gc_file` | path | no | – |
| `--analyse_TSDs` | flag | no | off |
| `--min_tsd` | integer | no | 3 |
| `--tsd_test_significance` | flag | no | off |
| `--analyse_palindromes` | flag | no | off |
| `--palindrome_min_arm` | integer | no | 3 |
| `--palindrome_max_spacer` | integer | no | 2 |
| `--palindrome_max_mismatches` | integer | no | 0 |
| `--analyse_bendability` | flag | no | off |
| `--analyse_melting` | flag | no | off |

## Output
Everything is written to `--outdir`, with file names starting with `--mge_name`:

```
outdir/
├── mge_name_YYYYmmdd-HHMMSS.log
├── files/
└── graphs/
```

Elements with fewer than 100 unique left or right flanks are skipped. You can edit the variable `MIN_UNIQUE_FLANKS` in `flankalyser.py`, though it is not recommended to go below 100 as results may not be very trustworthy.

### Files
| File | Produced | Contents |
|---|---|---|
| `MGE_FLANK_COUNTS.tsv` | always | Number of flank pairs in the input, unique left and right flanks, and unique left/right pairs |
| `MGE_SIGNIFICANT_POSITIONS.tsv` | always | Per position: p-value for the base composition differing from expected, whether it is significant after Bonferroni correction, and observed minus expected proportion of each base |
| `MGE_LEFT_DINUC_ENRICHMENT.tsv`, `MGE_RIGHT_DINUC_ENRICHMENT.tsv` | always | 4×4 table of log2 enrichment of each dinucleotide in the 10 bp nearest the insertion site vs. the 10 bp farthest away |
| `MGE_TSD_SUMMARY.txt` | `--analyse_TSDs` | Number and percentage of unique insertion sites with a TSD, breakdown by TSD length, and (with `--tsd_test_significance`) enrichment and p-value against shifted-boundary background |
| `MGE_PALINDROME_ENRICHMENT_CURVE_arm{A}_mm{M}_sp{S}.tsv` | `--analyse_palindromes` | Per minimum arm length: real vs. shuffled count of sequences with an inverted repeat at least that long, and the enrichment ratio. `A`, `M` and `S` are the min arm length, max mismatches and max spacer used |

### Graphs
| Graph | Produced | Shows |
|---|---|---|
| `MGE_DNA_FREQ.png` | always | Proportion of each base at each position, with significant positions marked |
| `MGE_LEFT_LOGO.png`, `MGE_RIGHT_LOGO.png` | always | Sequence logo of each flank |
| `MGE_DINUC_ENRICH.png` | always | Heatmaps of the dinucleotide enrichment tables |
| `MGE_TSD_LENGTH_DIST.png` | `--analyse_TSDs --tsd_test_significance` | TSD length distribution, observed vs. shifted-boundary background |
| `MGE_PALINDROME_COUNTS_arm{A}_mm{M}_sp{S}.png` | `--analyse_palindromes` | Real vs. shuffled counts by minimum arm length, with the inflection point marked |
| `MGE_BEND_PROFILE.png` | `--analyse_bendability` | Mean local bend per position, real vs. dinucleotide-shuffled |
| `MGE_MELTING_PROFILE.png` | `--analyse_melting` | Mean stacking ΔG per position, real vs. dinucleotide-shuffled |

## Examples
