## About
Analyse + visualise the regions surrounding insertion sites
## How to use
### Prerequisites
Needs output from atb - conversion file from lexicmap to flankalyser-compatible integrated (or just rewrite?)
### Set-up
conda env - mostly biopython stuff + graphing libs like matplotlib
### Options


### What I've changed
- I'm removing the GC% and species stuff for now - these should be part of it as they are doable from just the atb data, but they aren't directly related to the core functionality - make separate folders for these, including setup and how exactly to extract that information from atb, and how to pipe that output back into flankalyser (which would be - species/genus level flank data as well as more accurate deviations, instead of just deviations from 25%)
- the old (not great) palindrome analysis - while I do like the simplicity, the new find_pals.py is superior because you can adjust palindrome length, spacers, mismatches etc.
- same for tsd_detection.py vs tsd_analysis.py - they seem like they do the same thing but I have not used tsd_analysis.py in a bit and it seemed much slower?
- nucl_enrichment needs some work I think - should be able to pass in length of comparison? 25 seems to lenient, generally this should be TSD + a little bit of buffer
