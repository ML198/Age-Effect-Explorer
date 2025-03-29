library(tidyverse)
here::i_am("scripts/update_metadata.R")
data_dir <- "data"
raw_data_dir <- "data/raw_data/"
# metadata_path <- file.path(raw_data_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")
metadata_path <- file.path(raw_data_dir, "gtex_v10_metadata_exact.csv")

if (!file.exists(metadata_path)) stop("Error: Metadata file missing at ", metadata_path)

metadata <- read.csv(metadata_path)

# metadata <- metadata %>%
#   rename(donor = SUBJID, sex = SEX, age = AGE, death_type = DTHHRDY) %>%
#   mutate(
#     age_plot = as.numeric(sub("-.*", "", age)),
#     sex_plot = ifelse(sex == 1, "Male", "Female"))
# metadata$sex <- factor(metadata$sex)

save_path <- file.path(data_dir, "Updated_gtex_v10_metadata_exact.csv")
write.table(metadata, save_path, sep = "\t", row.names = FALSE)

