Вот текст в raw формате для прямого копирования в поле Edit:

# Genetic Analysis Scripts

This repository contains a collection of scripts for comprehensive genetic analysis and serves as supplementary material for our research paper, which is available as a preprint on medRxiv: https://doi.org/10.1101/2023.11.20.23298462

## Repository Contents

The repository includes the following R scripts for genetic analysis:

### `matching.R`
Script for performing case-control matching, including:
- Global population PCA-based matching
- Population-specific matching (European, African, Ad Mixed American)
- Quality assessment of matching using built-in metrics from the MatchIt package
- Visualization of matching results

### `functions.R`
Core functions library containing implementations of:
- HWE (Hardy-Weinberg Equilibrium) calculations
- CVAS (Common Variant Association Study) analysis functions
- RVAS (Rare Variant Association Study) analysis tools
- Statistical tests (Fisher's exact test, linear regression)
- Various plotting functions for result visualization

### `assoc.R`
Association analysis script implementing:
- CVAS analysis for missense variants and LoF mutations
- RVAS analysis with multiple testing approaches
- QQ-plot generation for case-control matching consistence.
- Meta-analysis utilities

## Analysis Components

### Clustering Analysis
Scripts for genetic data clustering that enable population structure identification and subgroup detection based on genetic markers.

### Case-Control Matching
Tools for matching control samples to study cases while accounting for various parameters, including demographic and genetic characteristics.

### Association Analysis
- **CVAS (Common Variant Association Study)**: Scripts for analyzing associations of common genetic variants
- **RVAS (Rare Variant Association Study)**: Tools for studying associations of rare genetic variants

## Related Publications

Detailed methodology and research findings are available in our further paper:

**Preprint**: https://doi.org/10.1101/2023.11.20.23298462 (prioritise the citation of the main article over the preprint)

## Citation

When using these scripts in your research, please prefer citation of our main paper.
