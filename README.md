Age-Effect-Explorer
================
2025-08-16

- [Age-Effect-Explorer](#age-effect-explorer)
- [Get start](#get-start)
- [Data Sources](#data-sources)
- [Project Structure](#project-structure)

## Age-Effect-Explorer

A Shiny application to visualize age- and sex-associated changes in gene
expression across 54 human tissues using GTEx v10 data.

## Get start

You can access the deployed Shiny app here:
<https://menghui.shinyapps.io/Age_Effect_Explorer/>

You can also reproduce the shiny app with
<https://github.com/menghui-c/Age-Effect-Explorer-Dashboard>, there has
instructions to install the app in details.

## Data Sources

- GTEx v10 TPM matrices, metadata (see [GTEx
  Portal](https://gtexportal.org/home/)).

## Project Structure

- `app.R`: Main Shiny application (UI + server logic).

- `scripts/p_value.R`: Runs linear regression analysis per tissue
  (log2TPM ~ age/sex/covariates), producing p-value tables.

- `scripts/GSEA.Rmd`: Performs Gene Set Enrichment Analysis (gsea) on
  ranked gene lists by tissue, stratified by age effect direction.

- `scripts/tissue_gene_extraction.R`:Extracts the list of available
  tissues and genes from GTEx data; writes tissue_names.txt and
  gene_names.txt.

- `scripts/update_metadata.R`:Cleans and standardizes the raw GTEx
  metadata.

- `scripts/covariated_download.R`:Downloads and unpacks GTEx covariates
  for regression modeling.

- `scripts/pval_distribution.R`:Aggregates p-value results across
  tissues and plots p-value distribution histograms

- `age_sig_genes_by_tissue.py`: Produces a summarized table including
  all the information of age-related significant genes across all 54
  tissues, each sheet represents one tissue.
