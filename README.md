# Schizophrenia glial program analysis

This repository contains the BA46 discovery analysis, BA11 validation, and an analysis of public bulk RNA-seq data from microglia-conditioned medium (MCM) experiments.

| Analysis | Notebooks | Data |
| --- | --- | --- |
| BA46 discovery: astrocyte–microglia module correlations, SNAP-a and GWAS | [primary_BA46](notebooks/primary_BA46/) | [primary_cohort](data/primary_cohort/) |
| BA11 validation: GSE254569 module transfer and de novo WGCNA | [validation_BA11](notebooks/validation_BA11/) | [validation_cohort/GSE254569](data/validation_cohort/GSE254569/) |
| MCM perturbation: mouse and human astrocyte scores | [perturbation](notebooks/perturbation/) | [perturbation_validation](data/perturbation_validation/) |

[Figure index](FIGURES.md) lists the numbered figures and their original names. [Data inventory](data/README.md) describes the inputs, intermediate data, and result tables. Reusable code lives in [sources](sources/), and WGCNA R scripts are in [Rcodes](Rcodes/).

The MCM review notebook shows the two astrocyte datasets. See the [perturbation README](notebooks/perturbation/README.md) for the commands to rebuild figures 24–25 or recompute the program scores.

Licensed under the [MIT License](LICENSE).
