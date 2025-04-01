library(dplyr)
library(tidyverse)
library(here)
library(broom)
library(qvalue)
# Load tissue names and metadata
app_dir <- here::here()
data_dir <- file.path(app_dir, "data")
covariates_dir <- file.path(data_dir, "covariates")
p_value_dir <- file.path(data_dir, "p_value")


tissue_file <- file.path(data_dir, "tissue_names.txt")
tissues <- readLines(tissue_file)

metadata.path <- file.path(data_dir, "Updated_gtex_v10_metadata_exact.csv")
metadata <- read.table(metadata.path, sep = "\t", header = TRUE) %>% select(donor, sex, age)

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
  
  # Extract sample columns and pivot to long format
  sample_cols <- colnames(exp)[-(1:2)]
  exp_long <- exp %>%
    select(Name, Description, all_of(sample_cols)) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "TPM") %>% 
    # Add donor information from metadata
    mutate(donor = str_replace_all(str_extract(sample, "^[^.]+\\.[^.]+"), "\\.", "-")) %>% 
    select(donor, Description, TPM) %>% 
    # merge with metadata
    left_join(metadata, by = "donor") 
  
  covariates.path <- file.path(covariates_dir, sprintf("%s.v10.covariates.csv", gsub(" ", "_", gsub("\\b([a-z])", "\\U\\1", tissue, perl = TRUE))))
  
  pval_df <- exp_long %>%
    group_by(Description) %>%
    summarise(
      p_value = {
        df_sub <- cur_data() %>%
          mutate(log2TPM = log2(TPM + 1), sex = as.factor(sex))
        
        if (file.exists(covariates.path)) {
          covariates <- read.table(covariates.path, sep = "\t", header = TRUE)
          df_sub <- left_join(df_sub, covariates, by = "donor")
        }
        
        if (sd(df_sub$log2TPM) < 1e-6) {
          NA_real_
        } else {
          remaining_vars <- setdiff(names(df_sub), c("donor", "Description", "TPM", "log2TPM"))
          
          if (length(remaining_vars) == 0 || nrow(df_sub %>% na.omit()) == 0) {
            NA_real_
          } else {
            fit <- tryCatch({
              lm(as.formula(paste("log2TPM ~", paste(remaining_vars, collapse = " + "))),
                 data = df_sub %>% na.omit()
              )
              # lm(log2TPM ~ age + sex, data = df_sub %>% na.omit())
            }, error = function(e) {
              print(paste("Error for gene:", first(df_sub$Description), "- Error Message:", e$message))
              NULL
            })
            
            if (is.null(fit)) {
              NA_real_
            } else {
              tidy_res <- tryCatch({
                tidy(fit)
              }, error = function(e){
                print(paste("Error during tidy for gene:", first(df_sub$Description), "- Error Message:", e$message))
                NULL
              })
              
              if(is.null(tidy_res)){
                NA_real_
              } else{
                age_pval <- tidy_res %>% filter(term == "age") %>% pull(p.value)
                if (length(age_pval) == 0) NA_real_ else age_pval
              }
            }
          }
        }
      },
      .groups = "drop"
    ) %>%
    rename(Gene = Description) %>%
    arrange(p_value) %>%
    mutate(
      Rank = row_number(),
      BH_adjusted = formatC(signif(p.adjust(p_value, method = 'BH'), 4), format = "e", digits = 4),
      `Storey's q-value` = formatC(signif(qvalue(p_value)$qvalues, 4), format = "e", digits = 4),
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

