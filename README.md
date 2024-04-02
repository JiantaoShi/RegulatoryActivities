---
title: "Scripts for exploring regulatory potential of predefined regions"
output:
  html_document:
    toc: true
---
***

# making windows regarding gene structures
```
usage: makeWindows.py [-h] -A ANNOTATION_FILE [-HWP WINDOW_PROMOTER] [-HWE WINDOW_ENHANCER] [-BN BIN_NUMBERS] [-C CUT_LENGTH] [-L MIN_LENGTH] -O OUTPUT_TAG

optional arguments:
  -h, --help            show this help message and exit
  -A ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        Gene annotation file.
  -HWP WINDOW_PROMOTER, --window_promoter WINDOW_PROMOTER
                        Half window size for promoters [100].
  -HWE WINDOW_ENHANCER, --window_enhancer WINDOW_ENHANCER
                        Half window size for enhancers [10000].
  -BN BIN_NUMBERS, --bin_numbers BIN_NUMBERS
                        Number of bins for gene bodies [20].
  -C CUT_LENGTH, --cut_length CUT_LENGTH
                        Numober of bases to cut from gene bodies [1000].
  -L MIN_LENGTH, --min_length MIN_LENGTH
                        Minimal length of gene body [1000].
  -O OUTPUT_TAG, --output_tag OUTPUT_TAG
                        Prefix for output files.
```


