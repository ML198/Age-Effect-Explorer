library(dplyr)
library(tidyverse)
library(here)
library(broom)
library(qvalue)

# Load tissue names and metadata
app_dir <- here::here()
data_dir <- file.path(app_dir, "data")
p_value_dir <- file.path(data_dir, "p_value")


tissue_file <- file.path(data_dir, "tissue_names.txt")
tissues <- readLines(tissue_file)
metadata.path <- file.path(data_dir, "Updated_GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
metadata <- read.table(metadata.path, sep = "\t", header = TRUE)

# Function to calculate p-values for all genes in a given tissue
calc_all_genes_pvalue_for_tissue <- function(tissue) {
  
  # Dynamically generate the URL for the tissue
  tissue_url <- sprintf("https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_%s.gct.gz", gsub(" ", "_", tissue))
  # Try to read the data from the URL
  exp <- tryCatch({
    temp_file <- tempfile(fileext = ".gz")
    download.file(tissue_url, temp_file, mode = "wb", quiet = TRUE)  # Set quiet = TRUE
    read.table(gzfile(temp_file), sep = "\t", skip = 2, header = TRUE)
  }, error = function(e) {
    return(NULL)
  })
  
  # Check the result
  if (is.null(exp)) {
    print("Error: Unable to download or read the file.")
  } 
  
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
  
  pval_df <- exp_long %>%
    group_by(Description) %>%
    summarise(
      p_value = {
        df_sub <- pick(everything())
        df_sub$log2TPM <- log2(df_sub$TPM + 1)
        
        # Check for variability
        if (sd(df_sub$log2TPM) == 0) {
          NA_real_
        } else {
          # Choose the appropriate model
          fit <- if (length(unique(df_sub$sex_plot)) == 1) {
            lm(log2TPM ~ age_plot, data = df_sub)
          } else {
            lm(log2TPM ~ age_plot + sex_plot, data = df_sub)
          }
          
          # Extract p-value using broom::tidy
          tidy_res <- tidy(fit)
          age_pval <- tidy_res %>% filter(term == "age_plot") %>% pull(p.value)
          
          if (length(age_pval) == 0) NA_real_ else age_pval
        }
      },
      .groups = "drop"
    ) %>%
    rename(Gene = Description) %>%
    arrange(p_value) %>%
    mutate(
      Rank = row_number(),
      BH_adjusted = formatC(signif(p.adjust(p_value, method = 'BH'), 4), format = "e", digits = 4),
      # `Storey's q-value` = ifelse(requireNamespace("qvalue", quietly = TRUE), 
      #                             formatC(signif(qvalue(p_value)$qvalues, 4), format = "e", digits = 4), 
      #                             NA_real_),
      `Storey's q-value` = ifelse(requireNamespace("qvalue", quietly = TRUE), 
                                  formatC(signif(qvalue(pval_df$p_value)$qvalues, 4), format = "e", digits = 4), 
                                  NA_real_),
      
      p_value = formatC(signif(p_value, 4), format = "e", digits = 4)
    ) %>% 
    select(Rank, Gene, p_value, BH_adjusted, `Storey's q-value`)
  
  return(pval_df)
}

# Loop through tissues and process each one
for (tissue in tissues) {
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
  output_path <- file.path(here::here("data/p_value"), paste0(gsub(" ", "_", tissue), "_pvalue_results.csv"))
  
  tryCatch({
    write.csv(result, output_path, row.names = FALSE)
    cat("Saved p-value results for tissue:", tissue, "\n")
  }, error = function(e) {
    cat("Error saving results for", tissue, ":", e$message, "\n")
  })
}

