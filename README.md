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
All flanks must be of the same length, but may be padded with `N`. Every flank must have a partner, else that ID will be skipped. Note the `--flank-len` argument must be accurate for the actual length of the flanks, else all input will be skipped.

```tsv
# Valid input - N and invalid characters are ignored
accession_1    NNAGCT    ACCNNN
accession_2    NNNNNN    AGTTTT
accession_3    AGNCCC    TTAANA
accession_4    HELLO1    world2 # characters not in [A,G,C,T] are converted to N
accession_5    agcttt    tgcnnn # lower-case nucleotides are fine
```
```tsv
# Invalid input - will log a warning and the offending row will be ignored
accession_1    NNAGCT # incomplete pair
accession_2    ACG    AGTTTT # incorrect length
accession_3    AGNCCC    AGCTTA|ATTGGG # incomplete pair
```
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
## Output

## Examples
