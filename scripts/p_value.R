library(dplyr)
library(tidyr)
library(here)

tissues <- readLines(here::here("gtex_v10_shiny/data/tissue_names.txt"))

metadata.path <- here::here("gtex_v10_shiny/data/raw_data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
metadata <- read.table(metadata.path, sep="\t", header=TRUE)
colnames(metadata) <- c("donor","sex","age","death_type")
metadata$age_plot <- sapply(metadata$age, function(a) as.numeric(strsplit(a, "-")[[1]][1]))
metadata$sex_plot <- ifelse(metadata$sex==1, "Male", "Female")

calc_all_genes_pvalue_for_tissue <- function(tissue) {
  exp.path <- sprintf(here::here("gtex_v10_shiny/data/raw_data/gene_tpm_v10_%s.gct.gz"), 
                      gsub(" ","_", tissue))
  if(!file.exists(exp.path)) {
    warning("File not found: ", exp.path)
    return(NULL)
  }
  
  exp <- read.table(gzfile(exp.path), sep="\t", skip=2, header=TRUE, check.names=FALSE)
  colnames(exp) <- gsub("\\.", "-", colnames(exp))
  
  sample_cols <- colnames(exp)[-(1:2)]
  
  exp_long <- exp %>%
    select(Name, Description, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "TPM"
    )
  
  exp_long$donor <- sapply(strsplit(exp_long$sample, "-"), function(v) paste(v[1], v[2], sep="-"))
  
  exp_long <- exp_long %>% left_join(metadata, by="donor")
  
  exp_long <- exp_long %>% filter(!is.na(TPM), !is.na(age_plot))
  
  pval_df <- exp_long %>%
    group_by(Description) %>%
    summarise(
      p_value = {
        df_sub <- pick(cur_data())  # Update deprecated function with pick()
        df_sub$logTPM <- log(df_sub$TPM + 1)
        fit <- lm(logTPM ~ age_plot + sex_plot, data=df_sub)
        s <- summary(fit)$coefficients
        if("age_plot" %in% rownames(s)) {
          s["age_plot","Pr(>|t|)"]
        } else {
          NA_real_
        }
      }
    ) %>%
    rename(Gene = Description) %>%
    ungroup()
  
  pval_df <- pval_df %>%
    arrange(p_value) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, Gene, p_value)
  
  return(pval_df)
}


# Run the function for all tissues
for (tissue in tissues) {
  cat("Processing tissue:", tissue, "\n")
  
  # Calculate p-values
  result <- tryCatch({
    calc_all_genes_pvalue_for_tissue(tissue)
  }, error = function(e) {
    cat("Error processing", tissue, ":", e$message, "\n")
    return(NULL)
  })
  
  # Check if result is NULL
  if (is.null(result)) {
    next
  }
  
  # Define output path
  output_path <- file.path(here::here("gtex_v10_shiny/data"), paste0(gsub(" ", "_", tissue), "_pvalue_results.csv"))
  
  # Save the results as a CSV
  tryCatch({
    write.csv(result, output_path, row.names = FALSE)
    cat("Saved p-value results for tissue:", tissue, "\n")
  }, error = function(e) {
    cat("Error saving results for", tissue, ":", e$message, "\n")
  })
}
