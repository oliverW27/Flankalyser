<p align="center">
  <img src="logo.png" width="180" alt="flankalyser logo">
</p>
<h1 align="center">Flankalyser</h1>
<p align="center">Analysis of the DNA flanking bacterial transposable elements</p>


## Quick start
Flankalyser looks for patterns around the insertion site of a (mobile) genetic element, using its flanking sequences from many genomes. By default it reports the base composition at each position (with significance), a sequence logo and dinucleotide enrichment near the insertion site. Optionally, it can also detect target site duplications (TSDs), look for palindromes, and profile DNA bendability and melting.

You need:
- the conda environment (see [Dependencies](#dependencies))
- a tab-separated file of left and right flanks for one element (see [Input](#input))
- This codebase

```bash
python flankalyser.py --mge_name IS26 \
  --summary_file examples/real_examples/flanks_IS26_subset.txt \
  --flank_len 50 --outdir results \
  --analyse_TSDs --analyse_palindromes \
  --analyse_bendability --analyse_melting
```

Graphs are written to `results/graphs/`, tables to `results/files/`, and a log of the run to `results/`. See [Output](#output) for details.

To run every `.txt` flank file in a folder, use `run_flankalyser_multi.sh` from the repository root, giving the folder, the flank length and the output directory:
```bash
bash run_flankalyser_multi.sh examples/real_examples 50 results
```
Each file's name without `.txt` is used as its `--mge_name`, e.g. `flanks_IS26_subset.txt` gives outputs named `flanks_IS26_subset_...`.

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
The input file must contain the flanking contexts of a given MGE in the following, tab-delimited format:

```tsv
[unique ID]   [left flank 1]|[optional: left flank 2]|...   [right flank 1]|[optional: right flank 2]...
```
All flanks must be of the same length, but may be padded with `N`. Note the `--flank_len` argument must be accurate for the actual length of the flanks, else all input will be skipped.

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

**Problematic input** (logs a warning):
```tsv
accession_1    NNAGCT
accession_2    ACG    AGTTTT
accession_3    AGNCCC    AGCTTA|ATTGGG
```
- `accession_1`: only one flank column, so the whole row is skipped.
- `accession_2`: the left flank is the wrong length, so it is skipped. The right flank is still used in the analyses that look at each flank on its own (nucleotide composition, sequence logo, dinucleotide enrichment, bendability, melting), but not in the TSD or palindrome analyses
- `accession_3`: different numbers of left and right flanks. All three flanks are used in the analyses that look at each flank on its own. The TSD and palindrome analyses pair the flanks in order (`AGNCCC` with `AGCTTA`) and leave out the unpaired `ATTGGG`.

The expected way of data generation for this tool is using [LexicMap](https://github.com/shenwei356/LexicMap) to query a database of assemblies, such as [AllTheBacteria](https://allthebacteria.org/).
A script `convert_lexicmap_fmt.py` is provided in `helpers/` to convert the output format to one compatible with Flankalyser.

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

### GC background file (`--gc_file`, optional)
The nucleotide composition test (step 2 of the [analysis](#analysis-steps)) asks whether the base proportions at each position differ from what is *expected*. By default, Flankalyser expects 25% of each base. Most bacterial genomes are not 50% GC, so in a genome with 35% GC every position will look AT-rich and come out significant, even with no target-site preference at all (see the `flanks_low_AT` and `flanks_high_AT` [toy examples](#toy-examples)). The GC file gives the expected GC content of the genomes the element is found in. The expected proportions then become A = T = (100 − GC)/2 % and C = G = GC/2 %, both for the significance test and for the baseline in the `DNA_FREQ` graph.

The file is tab-separated with a header line. Flankalyser looks for the row whose **second column exactly matches `--mge_name`** and reads the GC % (0–100) from the **third column**. Any other columns are ignored. Elements not found in the file fall back to 25% per base, with a warning.

```tsv
family	element	mean_gc	std_gc	n_accessions
IS21	IS100	50.8	1.07	65531
IS3	IS1397	50.58	0.74	42529
...
```

`examples/real_examples/element_gc.tsv` has the values for the four [real examples](#real-examples).

**How the background was calculated** for the examples:
1. For each genome (assembly) that contains the element, calculate the GC % of the whole assembly, with all contigs joined and only A/C/G/T counted.
2. Take the mean of these values over all genomes carrying the element, counting each genome once, however many copies of the element it has.

For example, step 1 for one gzipped assembly:
```python
import gzip

def genome_gc(fasta_gz):
    gc = acgt = 0
    with gzip.open(fasta_gz, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                seq = line.strip().upper()
                gc += seq.count("G") + seq.count("C")
                acgt += sum(seq.count(b) for b in "ACGT")
    return 100 * gc / acgt
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

## Analysis steps
Steps 1-3 always run; steps 4-7 are switched on with their `--analyse_xxx` option. In all steps, only **unique flanks** are counted even if they appear in multiple rows under different IDs. This is done as not every context represents a unique insertion event (could be due to e.g., resequencing of the same isolate, or clonal expansion). The exception to this is for TSD and palindromic analysis, here **unique left/right pairs** are used.
Several steps compare the real flanks with **dinucleotide-shuffled** copies of them. Shuffling destroys the order of the sequence but keeps its exact dinucleotide counts (after Altschul & Erickson, 1985), so the shuffled copies have the same base composition as the real flanks.

### 1. Input parsing
- Bases are converted to upper case, and any character other than A/C/G/T/N is replaced with N.
- Flanks that aren't `--flank_len` long are skipped with a warning.
- Elements with fewer than 100 unique left or right flanks are skipped (see [Output](#output)).
- With `--exclude_n`, flanks containing N are left out of steps 2 and 3. The other steps always leave them out de facto.

### 2. Nucleotide composition
- Counts A, C, G and T at each position of the left and right flanks (Ns are not counted).
- Tests each position with a chi-square goodness-of-fit test against the expected base proportions: 25% each by default, or A = T = (1 − GC)/2 and C = G = GC/2 if the element's GC content is given with `--gc_file`.
- Significance is Bonferroni-corrected for the number of positions tested (2 × `--flank_len`). In the graph, `*`, `**` and `***` mark corrected p < 0.05, 0.01 and 0.001.

### 3. Sequence logo and dinucleotide enrichment
- **Sequence logo:** the base frequencies from step 2, drawn with [logomaker](https://logomaker.readthedocs.io/).
- **Dinucleotide enrichment:** for each flank, compares how often each of the 16 possible dinucleotides occurs in the 10 bp nearest the insertion site vs. the 10 bp farthest from it. The result is log2(inner frequency / outer frequency), with a pseudocount of 0.5 per dinucleotide; positive values mean the dinucleotide is more common near the insertion site.

### 4. TSDs (`--analyse_TSDs`)
- For each unique left/right flank pair, finds the longest exact match between the end of the left flank and the start of the right flank, from `--min_tsd` up to 25 bp. A TSD can't contain N.
- **Significance (`--tsd_test_significance`):** the same search is repeated with the boundary moved 10–25 bp into the left flank and 10–25 bp into the right flank (32 shifted boundaries). Matches found there are not real TSDs, so the average number of sites with a match per shifted boundary is the background. The observed number is compared to it with a one-sided binomial test, and reported as an enrichment (observed ÷ background) and p-value.

### 5. Palindromes (`--analyse_palindromes`)
- Joins each unique left/right flank pair into one sequence (pairs containing N are skipped) and finds the longest palindrome in it. Every position is tried as the centre, with a spacer of 0 to `--palindrome_max_spacer` bp between the two arms, and arms are extended outwards while the bases are complementary, allowing up to `--palindrome_max_mismatches` mismatched pairs.
- For each arm length from `--palindrome_min_arm` to 50 bp, counts how many sequences have an inverted repeat at least that long, in the real sequences and in 100 dinucleotide-shuffled copies of them.
- **Enrichment** at each arm length is the real count ÷ the mean shuffled count. The **inflection point** is the shortest arm length with at least 2× enrichment, i.e. the length from which palindromes are clearly more common than chance. This metric is somewhat arbitrary and should only be interpreted together with the nucleotide distribution graph and palindrome distribution and count.

### 6. DNA bendability (`--analyse_bendability`)
- Uses the trinucleotide nucleosome-positioning model of Goodsell & Dickerson (1994). Each trinucleotide is given a roll angle (0–10°) according to how strongly it prefers to face away from the histone core, with a constant twist of 34.3° per base step.
- The base pair frames are built up step by step along each flank, and the local bend at a position is the angle between the base pair normals one step either side of it.
- The mean bend at each position is plotted for the real flanks and for one dinucleotide-shuffled copy of each flank. Flanks containing N are skipped, and at most 5000 randomly chosen flanks per side are used. There is no significance test.

### 7. DNA melting (`--analyse_melting`)
- Gives each dinucleotide step its nearest-neighbour stacking free energy (ΔG°37, SantaLucia 1998 unified parameters). Less negative ΔG means the DNA is easier to melt.
- The mean ΔG at each position is plotted for the real flanks (with a standard-error band) and for dinucleotide-shuffled copies. Flanks containing N are skipped.

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
### Toy examples
`helpers/make_synth_data.py` is provided to make synthetic data - five example inputs of various target site patterns are present, but playing around with different types of preference patterns may be helpful to understand real results. The number and length of the flanks can be changed with `N` and `FLANK_LEN` at the top of the script. The sequences are random, so exact numbers will differ slightly between runs.

```bash
for f in examples/toy_examples/*.txt; do
    python flankalyser.py --mge_name $(basename $f .txt) --summary_file $f --flank_len 50 \
        --outdir examples/toy_examples --analyse_TSDs --tsd_test_significance \
        --analyse_palindromes --analyse_bendability --analyse_melting
done
```

| Dataset | How it is made | What it shows |
|---|---|---|
| `flanks_no_preference` | Every base random (25% each) | At every position, each nucleotide is equally common - there is NO target site motif. This is a baseline against which to compare any other results!
| `flanks_low_AT`, `flanks_high_AT` | Every base drawn with 20% or 75% A+T | A composition bias across the whole flank. These datasets also have more chance TSDs (around 5% of sites vs. 1.5%), but `--tsd_test_significance` then shows that the background is just as high. |
| `flanks_strong_short_seq` | `CTAG` directly left of the insertion site in 80% of flanks | A strong, short target-site preference: positions −4 to −1 are highly significant, the left logo spells `CTAG`, and `CT`, `TA` and `AG` are about 2× enriched near the insertion site in the left flank only. |
| `flanks_weak_long_seq` | 15 bp motif `GCCGATCGCATCGGC` spanning the insertion site in 10% of flanks | A weak, long preference: most motif positions are still significant with 1,000 flanks, despite the motif being present in only 10% of them. The motif is also an inverted repeat (`GCCGAT`, 3 bp spacer `CGC`, `ATCGGC`), so it isn't found with the default `--palindrome_max_spacer 2`. With `--palindrome_max_spacer 3`, palindromes are about 2.3× enriched from an arm length of 6 bp. |

### Real examples
`examples/real_examples/` contains four real datasets, each a random subset of 1,000 unique insertion sites (left/right flank pairs, 50 bp) from flanks extracted from [AllTheBacteria](https://allthebacteria.org/) with the workflow above. Whole rows were kept, so genomes sharing the same insertion site are still in the files. The shuffled controls are random, so palindrome, bendability and melting results will differ slightly between runs.

```bash
for n in IS1397 IS100 IS26 IS91; do
    python flankalyser.py --mge_name $n \
      --summary_file examples/real_examples/flanks_${n}_subset.txt --flank_len 50 \
      --outdir examples/real_examples --analyse_TSDs --tsd_test_significance \
      --analyse_palindromes --analyse_bendability --analyse_melting \
      --gc_file examples/real_examples/element_gc.tsv 
done
```

| Dataset | What it shows |
|---|---|
| `flanks_IS1397_subset.txt` | A palindromic target ([a REP in E. coli](https://pmc.ncbi.nlm.nih.gov/articles/PMC169884/)). Long inverted repeats are strongly enriched. Short TSDs (mostly 3-4 bp) are found at 32% of sites, as has been [reported](https://pmc.ncbi.nlm.nih.gov/articles/PMC95333/). |
| `flanks_IS100_subset.txt` | A strong TSD with little target sequence preference. 5 bp TSDs (and some 6 bp) are found at 68% of sites. 10-bp AT peaks are visible in the bending graph, which indicates it inserts into bent DNA |
| `flanks_IS26_subset.txt` | An 8 bp TSD at only a small fraction of sites - IS26 rearranges its insertion contexts which disrupts these TSDs. Many positions near the site are significant, but no base is more than 15% above expected, so the target preference is weak. |
| `flanks_IS91_subset.txt` | A target motif without a TSD, this element does not create them. IS91 has a strong preference to [insert next to the GAAC/CAAG motif](https://pubmed.ncbi.nlm.nih.gov/2552258/), seen here at positions +1 to +4 in the right-flank sequence logo and nucleotide distribution graph |
