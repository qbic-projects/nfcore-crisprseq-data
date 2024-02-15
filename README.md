# nf-core/crisprseq data

This repository contains the scripts and data used in the nf-core/crisprseq publication.

## Benchmarking of resources

- `benchmarking_results/`:
    Contains tables of resources of different runs with nf-core/crisprseq and CRISPR-Analytics.

- `crisprseq_vs_crispra.R`:
    R script used to calculate the average run time and CPU hours of both workflows.

- `public_data/`:
    Contains the projects from the European Nucleotide Archive (ENA) repository filtered by the keyword `crispr`.

- `analysis_public_data.R`:
    R script used to classify public available projects into different categories, calculate the percentages and plot them.