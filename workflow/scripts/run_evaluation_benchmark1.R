#!/usr/bin/env Rscript --vanilla

# Benchmarking: evaluate performance across 4 test cases:
# 1. translation_strict : TP = translation,  TN = background+abundance+buffering
# 2. translation_only   : TP = translation,  TN = background only (abundance & buffering excluded)
# 3. any_regulated      : TP = translation+abundance+buffering, TN = background
# 4. buffering          : TP = buffering, TN = background+translation+abundance

# Input:
# identifier : gene ID
# class      : predicted class (Riborex - only translation)
# effect     : log2FC or equivalent effect size
# score      : -log10(adjpval)
