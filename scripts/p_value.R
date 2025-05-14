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
# tissues <- (c("cervix ectocervix", "cervix endocervix", "fallopian tube", "ovary", "prostate", "testis", "uterus", "vagina"))
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
    left_join(metadata, by = "donor") %>% 
    group_by(sex, Description)
    
  
  covariates.path <- file.path(covariates_dir, sprintf("%s.v10.covariates.csv", gsub(" ", "_", gsub("\\b([a-z])", "\\U\\1", tissue, perl = TRUE))))
  
  pval_df <- exp_long %>%
    mutate(log2TPM = log2(TPM + 1), sex = as.factor(sex)) %>%
    group_by(Description) %>%
    group_modify(~ {
      df_sub <- .x 
      
      # df_sub <- exp_long %>%
      #   filter(Description == "DDX3Y")  %>%
      #   mutate(log2TPM = log2(TPM + 1), sex = as.factor(sex))
      
      # merge covariates if available
      if (file.exists(covariates.path)) {
        covariates <- read.table(covariates.path, sep = "\t", header = TRUE)
        df_sub <- left_join(df_sub, covariates, by = "donor")
      }
      
      # skip if TPM is too constant
      if (nrow(df_sub) == 0 || sd(df_sub$log2TPM, na.rm = TRUE) < 1e-6) {
        return(tibble(p_value_age = NA_real_, p_value_sex = NA_real_, age_coef = NA_real_, sex_coef = NA_real_, sex_bias = NA_character_))
        # return(tibble(p_value = NA_real_, age_coef = NA_real_))
      }
      
      # Identify remaining variables for the regression model
      model_vars <- setdiff(names(df_sub), c("donor", "Description", "TPM", "log2TPM"))
      if (!all(c("age", "sex") %in% model_vars)) {
        return(tibble(p_value_age = NA_real_, p_value_sex = NA_real_, age_coef = NA_real_, sex_coef = NA_real_, sex_bias = NA_character_))
      }
      
      # Remove sex if it has only one level
      if (length(unique(df_sub$sex)) == 1) {
        model_vars <- setdiff(model_vars, "sex")
      }
      
      # If no remaining variables or if data has too many missing values, return NAs
      if (length(model_vars) == 0 || nrow(na.omit(df_sub)) == 0) {
        return(tibble(p_value_age = NA_real_, p_value_sex = NA_real_, age_coef = NA_real_, sex_coef = NA_real_, sex_bias = NA_character_))
      } 
      
      fit <- tryCatch({
        lm(as.formula(paste("log2TPM ~", paste(model_vars, collapse = " + "))), data = na.omit(df_sub))
      }, error = function(e) NULL)
      
      if (is.null(fit)) {
        return(tibble(p_value_age = NA_real_, p_value_sex = NA_real_, age_coef = NA_real_, sex_coef = NA_real_, sex_bias = NA_character_))
      } 
      
      # Extract tidy results from the model and handle errors
      tidy_res <- tryCatch({ tidy(fit) }, error = function(e) NULL)
      
      # If tidy results are missing, return NAs
      if (is.null(tidy_res)) {
        return(tibble(p_value_age = NA_real_, p_value_sex = NA_real_, age_coef = NA_real_, sex_coef = NA_real_, sex_bias = NA_character_))
      }
      p_value_age <- tidy_res %>% filter(term == "age") %>% pull(p.value)
      p_value_sex <- tidy_res %>% filter(term == "sexMale") %>% pull(p.value)
      age_coef <- coef(fit)["age"]
      age_sign <- sign(age_coef)
      sex_coef <- coef(fit)["sexMale"]
      sex_bias <- case_when(
        sex_coef > 0 ~ "Male",
        sex_coef < 0 ~ "Female",
        is.na(sex_coef) ~ NA_character_,
        TRUE ~ NA_character_
      )
      
      
      # Return p-value, coefficient for age, and sign of the age coefficient
      tibble(
        p_value_age = if (length(p_value_age) == 0) NA_real_ else p_value_age,
        p_value_sex = if (length(p_value_sex) == 0) NA_real_ else p_value_sex,
        age_coef = if (is.null(age_coef)) NA_real_ else age_coef,
        age_sign = if (is.null(age_sign)) NA_real_ else age_sign,
        sex_coef = if (is.null(sex_coef)) NA_real_ else sex_coef,
        sex_bias = sex_bias
      )
    },
    .groups = "drop"
    ) %>%
    ungroup() %>%
    rename(Gene = Description) %>%
    arrange(p_value_age) %>%
    mutate(
      Rank = row_number(),
      BH_adjusted_age = p.adjust(p_value_age, method = 'BH'),
      BH_adjusted_sex = p.adjust(p_value_sex, method = 'BH')
    ) %>%
    mutate(
      # `Storey's q-value` = formatC(signif(qvalue(p_value)$qvalues, 3), format = "e", digits = 3),
      p_value_age = formatC(signif(p_value_age,4), format = "e", digits = 4),
      BH_adjusted_age = formatC(signif(BH_adjusted_age, 4), format = "e", digits = 4),
      age_coef = formatC(signif(age_coef, 4), format = "e", digits = 4),
      
      p_value_sex = formatC(signif(p_value_sex,4), format = "e", digits = 4),
      BH_adjusted_sex = formatC(signif(BH_adjusted_sex, 4), format = "e", digits = 4),
      sex_coef = formatC(signif(sex_coef, 4), format = "e", digits = 4),
      
    ) %>%
    select(Rank, Gene, p_value_age, BH_adjusted_age, age_coef, age_sign, p_value_sex, BH_adjusted_sex, sex_coef, sex_bias)
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

