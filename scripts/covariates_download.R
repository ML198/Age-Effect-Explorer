library(tidyverse)
app_dir <- here::here()
data_dir <- file.path(app_dir, "data")
covariates_dir <- file.path(data_dir, "covariates")

url <- "https://storage.googleapis.com/adult-gtex/bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL_covariates.tar"
destfile <- file.path(covariates_dir, "GTEx_Analysis_v10_eQTL_covariates.tar")
download.file(url, destfile, mode = "wb")
cat("Download complete: ", destfile, "\n")

# Extract the contents
untar(destfile, exdir = covariates_dir)

# List the extracted files
list.files(covariates_dir)


# List all .txt files in the directory
txt_files <- list.files(covariates_dir, pattern = "\\.txt$", full.names = TRUE)

# Rename .txt to .csv
for (file in txt_files) {
  new_file <- sub("\\.txt$", ".csv", file)
  file.rename(file, new_file)
}
txt_files <- list.files(covariates_dir, pattern = "\\.csv$", full.names = TRUE)

cat("All .txt files have been renamed to .csv.")

# load the metadata
metadata.path <- file.path(data_dir, "Updated_gtex_v10_metadata_exact.csv")
metadata <- read.table(metadata.path, sep = "\t", header = TRUE) %>% select(donor, sex, age)

transpose_and_overwrite <- function(file) {
  data <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Transpose, add donor column, and reorder using dplyr
  transposed_data <- data %>%
    select(-1) %>% # Exclude the first column
    t() %>%
    as.data.frame() %>%
    mutate(donor = colnames(data)[-1]) %>%
    select(donor, everything())
  
  # Rename columns
  colnames(transposed_data) <- c("donor", as.character(data[, 1]))
  rownames(transposed_data) <- seq_len(nrow(transposed_data))
  
  # Get the overlapping column names (excluding "donor")
  overlapping_cols <- intersect(names(metadata), names(transposed_data))
  overlapping_cols <- overlapping_cols[overlapping_cols != "donor"]
  
  # Remove overlapping columns from covariates
  covariates_filtered <- transposed_data %>%
    select(-all_of(overlapping_cols))
  
  # Overwrite the file
  write.table(covariates_filtered, file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("File overwritten with transposed data:", file, "\n")
  dim(covariates_filtered)
}
# Apply the function to all .txt files
lapply(txt_files, transpose_and_overwrite)


# 
# # Extract tissue names from the covariate file names
# covariate_tissue_names <- gsub(".v10.covariates.csv", "", basename(txt_files))
# 
# # Clean up tissue names in 'tissues' to match the covariate file names
# tissues_cleaned <- tolower(gsub(" ", "_", tissues)) # Convert spaces to underscores and lowercase
# 
# # Check if each cleaned tissue in 'tissues_cleaned' has a corresponding covariate file
# missing_tissues <- setdiff(tissues_cleaned, tolower(covariate_tissue_names))
# 
# if (length(missing_tissues) == 0) {
#   print("All tissues have corresponding covariate files.")
# } else {
#   print("The following tissues are missing covariate files:")
#   print(missing_tissues)
# }
