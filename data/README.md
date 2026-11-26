# Data inventory

| Location | Contents |
| --- | --- |
| [primary_cohort/snRNA_donor](primary_cohort/snRNA_donor/) | BA46 donor-level single-nucleus expression |
| [primary_cohort/pb_ct](primary_cohort/pb_ct/) / [pb_ct_filtered](primary_cohort/pb_ct_filtered/) | BA46 cell type pseudobulk matrices before and after gene filtering |
| [primary_cohort/WGCNA_modules](primary_cohort/WGCNA_modules/) / [WGCNA_GO_enrichr](primary_cohort/WGCNA_GO_enrichr/) | BA46 module genes, preservation statistics and GO enrichment |
| [primary_cohort/corr_diff](primary_cohort/corr_diff/) / [gwas_overlap](primary_cohort/gwas_overlap/) | Module correlation differences and gene-level GWAS overlap |
| [validation_cohort/GSE254569](validation_cohort/GSE254569/) | BA11 input, per-donor and pseudobulk intermediates, WGCNA modules and GO enrichment |
| [SNAP_loadings](SNAP_loadings/) | Published human SNAP-a and SNAP-n loadings |
| [Homo_sapiens.GRCh38.114.gtf.gz](Homo_sapiens.GRCh38.114.gtf.gz) / [GWAS associations](gwas-association-downloaded_2026-03-23-MONDO_0005090.tsv) | Shared gene annotation and BA46 GWAS input |

The [perturbation_validation](perturbation_validation/) directory holds all MCM analysis data:

| Subdirectory | Contents |
| --- | --- |
| [GSE175485](perturbation_validation/GSE175485/) | Downloaded mouse bulk RNA-seq workbooks, including counts and log2FPKM |
| [GSE228927](perturbation_validation/GSE228927/) | Downloaded mouse astrocyte FPKM matrix for resting-, LPS- and hemin-MCM |
| [GSE282665](perturbation_validation/GSE282665/) | Downloaded human gene TPM matrices, split by recipient or experiment |
| [reference](perturbation_validation/reference/) | GEO sample metadata, human–mouse orthology and source hashes |
| [results](perturbation_validation/results/) | Sample manifest, program scores, mapping and coverage, condition summaries, paired differences and contrast statistics. Full tables include earlier exploratory comparisons. |

Several large cohort intermediates and shared references are available in the local workspace and excluded from Git. The numbered figure files are indexed in [FIGURES.md](../FIGURES.md).
