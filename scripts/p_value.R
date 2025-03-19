library(dplyr)
library(tidyverse)
library(here)

# Load tissue names and metadata
tissues <- readLines(here::here("gtex_v10_shiny/data/tissue_names.txt"))
metadata.path <- here::here("gtex_v10_shiny/data/raw_data/Updated_GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
metadata <- read.table(metadata.path, sep = "\t", header = TRUE)

# Function to calculate p-values for all genes in a given tissue
calc_all_genes_pvalue_for_tissue <- function(tissue) {
  exp.path <- sprintf(here::here("gtex_v10_shiny/data/raw_data/gene_tpm_v10_%s.gct.gz"), 
                      gsub(" ", "_", tissue))
  
  # Check if the file exists
  if (!file.exists(exp.path)) {
    warning("File not found: ", exp.path)
    return(NULL)
  }
  
  # Read the expression data
  exp <- read.table(gzfile(exp.path), sep = "\t", skip = 2, header = TRUE, check.names = FALSE)
  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  
  # Extract sample columns and pivot to long format
  sample_cols <- colnames(exp)[-(1:2)]
  exp_long <- exp %>%
    select(Name, Description, all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "TPM")
  
  # Add donor information from metadata
  exp_long$donor <- sapply(strsplit(exp_long$sample, "-"), function(v) paste(v[1], v[2], sep = "-"))
  exp_long <- exp_long %>% left_join(metadata, by = "donor")
  
  # Filter for valid TPM and age values
  exp_long <- exp_long %>% filter(!is.na(TPM), !is.na(age_plot))
  
  # Calculate p-values for each gene
  pval_df <- exp_long %>%
    group_by(Description) %>%
    summarise(
      p_value = {
        df_sub <- cur_data_all()
        df_sub$logTPM <- log(df_sub$TPM + 1)
        df_sub <- df_sub %>%
          mutate(sex_plot = factor(sex_plot),
                 age_plot = factor(age_plot))
        
        # Check for single level in sex_plot
        if (length(unique(df_sub$sex_plot)) == 1) {
          fit <- lm(logTPM ~ age_plot, data = df_sub)
        } else {
          fit <- lm(logTPM ~ age_plot + sex_plot, data = df_sub)
        }
        
        s <- summary(fit)$coefficients
        if("age_plot" %in% rownames(s)) {
          s["age_plot","Pr(>|t|)"]
        } else {
          NA_real_
        }
      }
    ) %>%
    rename(Gene = Description) %>%
    arrange(p_value) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Gene, p_value)
  
  return(pval_df)
}

# Loop through tissues and process each one
for (tissue in "cervix ectocervix") {
  cat("Processing tissue:", tissue, "\n")
  
  # Calculate p-values
  result <- tryCatch({
    calc_all_genes_pvalue_for_tissue(tissue)
  }, error = function(e) {
    cat("Error processing", tissue, ":", e$message, "\n")
    return(NULL)
  })
  
  # If no result, skip the tissue
  if (is.null(result)) {
    next
  }
  
  # Define output path and save results
  output_path <- file.path(here::here("gtex_v10_shiny/data"), paste0(gsub(" ", "_", tissue), "_pvalue_results.csv"))
  
  tryCatch({
    write.csv(result, output_path, row.names = FALSE)
    cat("Saved p-value results for tissue:", tissue, "\n")
  }, error = function(e) {
    cat("Error saving results for", tissue, ":", e$message, "\n")
  })
}
