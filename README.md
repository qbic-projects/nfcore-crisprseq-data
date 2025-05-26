# nf-core/crisprseq data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15519430.svg)](https://doi.org/10.5281/zenodo.15519430)

This repository contains the scripts and data used in the nf-core/crisprseq publication.

## Benchmarking of resources

- `benchmarking_results/`:
    Contains tables of resources of different runs with nf-core/crisprseq (`crisprseq_full`) and CRISPR-Analytics (`crispra_full`).
    Contains tables of analysis results run with nf-core/crisprseq (`crisprseq_spikes`). The tables are obtained from the MultiQC process, which summarises the results. Samples from two projects were analysed, Connelly JP et.al. (2019) (crisp.py) and Sanvicente-García M et.al. (2023) (crispr-A).

- `crisprseq_vs_crispra.R`:
    R script used to calculate the average run time and CPU hours of both workflows.

- `spikein_benchmarking.R`:
    R script used to plot a comparison of indels detected in spike-in samples with nf-core/crisprseq and compared with the original publications which analysed the data with crisp.py and CRISPR-A.

- `public_data/`:
    Contains the projects from the European Nucleotide Archive (ENA) repository filtered by the keyword `crispr`.

- `analysis_public_data.R`:
    R script used to classify public available projects into different categories, calculate the percentages and plot them.