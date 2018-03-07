# Data files

## input_data/ 

Files that do not get edited. Includes:
- NCI MAF file (tarred and gzipped), in NCI directory
- virusscan counts file (`vscan_counts.tsv`), for determining HPV status

## output_data/ 

Analysis output files. Includes:
- MAFs after split, removing DNP, changing to hg19, etc
- trinucleotide context file from the main MAF with no HPV split.

## Trinucleotide context file

You can load() this into R. It is a list structure, where:
- trinuc.contexts[[1]] gives you the summary used for the selection analysis (average of contexts)
- trinuc.contexts[[2]][[1]] gives you the output from the first tumor
- trinuc.contexts[[2]][[2]] gives the second, etc.