here::i_am("scripts/sig_gene.R")
library(tidyverse)
library(here)

# Directories
data_dir <- here::here("data")
pval_dir <- file.path(data_dir, "p_value")
out_xlsx <- file.path(pval_dir, "overlap_sig_age_sex.xlsx")
files <- list.files(pval_dir, pattern = "\\.csv$", full.names = TRUE)
files <- list.files(
  path = "data/p_value", 
  pattern = "_pvalue_results\\.csv$", 
  full.names = TRUE
)

# Function: count significant genes per tissue
count_sig <- function(fp) {
  df <- read_csv(fp, show_col_types = FALSE)
  
  if (!all(c("Gene", "BH_adjusted_age", "BH_adjusted_sex") %in% names(df))) {
    message("Skipping file (missing required columns): ", basename(fp))
    return(NULL)
  }
  
  tissue <- tools::file_path_sans_ext(basename(fp)) |>
    sub("_pvalue_results$", "", x = _) |>
    str_replace_all("_", " ") |>
    str_trim()
  
  tibble(
    Tissue = tissue,
    n_age_sig = sum(df$BH_adjusted_age < 0.05, na.rm = TRUE),
    n_sex_sig = sum(df$BH_adjusted_sex < 0.05, na.rm = TRUE)
  )
}
sig_counts <- map_dfr(files , count_sig)%>%
  mutate(
    label_tissue = if_else(n_sex_sig > 600 | n_age_sig > 1000, Tissue, NA_character_)
  )

# Plot: scatter of age vs sex significant genes per tissue
p <- sig_counts %>%
  ggplot(aes(x = log10(n_age_sig), y = log10(n_sex_sig), label = Tissue)) +
  geom_point(size = 3, color = "#2E86AB") +
  ggrepel::geom_text_repel(
    aes(label = label_tissue),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    segment.color = "grey50"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Significant Genes per Tissue",
    subtitle = "BH-adjusted p < 0.05 for Age and Sex effects",
    x = "Number of Age-significant genes (log10 scale)",
    y = "Number of Sex-significant genes (log10 scale)"
    ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13),
    panel.grid.minor = element_blank()
  )
p

find_overlap <- function(fp) {
  df <- read_csv(fp, show_col_types = FALSE)
  
  overlap <- df %>%
    filter(BH_adjusted_age < 0.05, BH_adjusted_sex < 0.05) 

  if (nrow(overlap) == 0) return(NULL)
  
  tissue <- tools::file_path_sans_ext(basename(fp)) |>
    str_remove("^GSEA[_-]?") |>
    str_replace_all("[_.]+", " ") |>
    str_trim()
  
  tibble(Tissue = tissue, overlap)
}
overlap_all <- map_dfr(files, find_overlap)


# Summary counts
counts <- overlap_all %>%
  count(Tissue, name = "n_overlap") %>%
  arrange(desc(n_overlap))

wb <- createWorkbook()

# summary sheet first
addWorksheet(wb, "Summary")
writeData(wb, "Summary", counts)

# tissue sheets
for (t in unique(overlap_all$Tissue)) {
  df_tissue <- overlap_all %>% filter(Tissue == t)
  addWorksheet(wb, str_sub(t, 1, 31))  # Excel sheet name limit
  writeData(wb, str_sub(t, 1, 31), df_tissue)
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

