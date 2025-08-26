
### Gene Set Enrichment ANalysis (GSEA)

#   This script runs GSEA on the age-associated genes in the GTEx 
# Consortium v10 dataset.
here::i_am("scripts/GSEA.Rmd")
knitr::opts_chunk$set(echo = TRUE)

suppressMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(patchwork)
})

# ---------- ggplot theme ----------
theme_ms <- function(base_size = 14, base_family = "") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1),
      strip.background = ggplot2::element_rect(linetype = "blank"),
      plot.title = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
      panel.grid.minor = ggplot2::element_line(colour = "grey90", linewidth = 0.5),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.5),
      complete = FALSE
    )
}

# ---------- paths & params ----------
app_dir    <- here::here()
data_dir   <- file.path(app_dir, "data")
pval_dir   <- file.path("/Users/mingruili/Desktop/untitled folder/csv")
plot_dir   <- file.path(app_dir, "figures")
out_dir    <- file.path(data_dir, "gsea")
tissue_file <- file.path(data_dir, "tissue_names.txt")
tissues <- readLines(tissue_file)
padj_thresh      <- 5e-2
top_n            <- 8

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)


#  ------------------- MSigDB genes (C2:CP:WIKIPATHWAYS)  -------------------
geneset <- msigdbr(
  species     = "Homo sapiens",
  category    = "C2",
  subcategory = "CP:WIKIPATHWAYS"
)

gene_list <- split(geneset$gene_symbol, geneset$gs_name)
pathway_map <- geneset %>%
  distinct(gs_name, gs_id, gs_description, .keep_all = FALSE) %>%
  transmute(
    pathway_key  = gs_name,
    pathway_id   = gs_id,
    pathway_desc = coalesce(gs_description,
                            str_replace_all(gs_name, "_", " "))
  )

#  ------------------- Main GSEA function  -------------------
# Define GSEA function to use `fgsea()` on predefined gene set

doGSEA <- function(pval_table, gene_list, tissue_name = "",
                   restrict_ranking = FALSE,
                   top_n = 8) {
  

  df <- pval_table %>%
    select(Gene, p_value_age, age_coef) %>%
    mutate(
      p_value_age = as.numeric(p_value_age),
      age_coef    = as.numeric(age_coef),
      rank        = sign(age_coef) * -log10(p_value_age)
    ) %>%
    filter(is.finite(rank), is.finite(p_value_age), !is.na(Gene)) %>% 
    arrange(desc(rank)) %>%            
    distinct(Gene, .keep_all = TRUE)
  
  set.seed(1234)
  Z <- setNames(df$rank, df$Gene)
  
  gsea_res <- fgseaMultilevel(pathways = gene_list, stats = Z, 
                              maxSize = 300, nPermSimple = 10000, nproc = 1) %>%
    arrange(padj) %>%
    mutate(pathway_key  = pathway,  
           Direction = ifelse(NES > 0, "Upregulated", "Downregulated")) %>% 
    filter(!is.na(NES)) %>% 
    dplyr::left_join(pathway_map, by = "pathway_key") %>%
    dplyr::mutate(plot_label = stringr::str_wrap(pathway_desc, width = 40))
  
  
  safe_blank_plot <- function(title = "") {
    ggplot2::ggplot() + ggplot2::theme_void() +
      ggplot2::labs(title = stringr::str_to_title(gsub("_"," ", title)))
  }
  
  if (nrow(gsea_res) == 0) {
    message(sprintf("[%s] gsea results were empty。", tissue_name))
    return(list(result = gsea_res, plot = safe_blank_plot(tissue_name)))
  }
  

  to_plot <- gsea_res %>%
    # mutate(pathway = str_wrap(gsub("_"," ", pathway), width = 40)) %>%
    group_by(Direction) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    droplevels()
  
  # compress to one strip if only one side requested
  p <- ggplot(to_plot, aes(x = -log10(padj),
                           y = reorder(plot_label, -padj),
                           size = size,
                           colour = abs(NES))) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               colour = "darkorange", linewidth = 0.8) +
    facet_grid(rows = vars(Direction), space = "free", scales = "free", drop = TRUE) +
    theme_ms() +
    labs(x = "-log10 FDR", y = "Pathway",
         title = stringr::str_to_title(tissue_name),
         colour = "|NES|", size = "Gene set size")
  
  list(result = gsea_res, plot = p)
}

# -------------------GSEA tables generating -------------------

make_and_save_tables <- function(gsea_res, tissue_name,
                                 padj_thresh = 5e-2,
                                 out_dir = "data/gsea") {
  
  outfile <- file.path(out_dir, sprintf("GSEA_%s.csv", tissue_name))
  
  if (!all(c("pathway_id", "pathway_desc") %in% names(gsea_res))) {
    gsea_res <- gsea_res %>%
      dplyr::left_join(pathway_map %>% dplyr::select(pathway_key, pathway_id, pathway_desc),
                       by = "pathway_key")
  }
  
  gsea_table_skeleton <- tibble::tibble(pathway_id = character(),pathway = character(), Direction  = character(),NES = double(),pval = double(),padj = double())
  
  if (nrow(gsea_res) == 0) {
    readr::write_csv(gsea_table_skeleton, outfile)  
    return(list(table = gsea_table_skeleton, file = outfile))
  }
  
  keep_df <- gsea_res %>%
    filter(padj < padj_thresh) %>%
    # mutate(pathway = gsub("_"," ", pathway)) %>%
    select(pathway_id, pathway = pathway_desc, Direction, NES, pval, padj)
  
  readr::write_csv(keep_df, outfile)
  
  list(table = keep_df, file = outfile)
}

# ------------------- Run across tissues -------------------
all_tables <- list()

for (t in tissues) {
  tissue_name <- gsub(" ","_",t)
  pval_table <- read.csv(here::here(pval_dir, sprintf("%s_pvalue_results.csv", tissue_name)), header = TRUE,   na.strings = c("", "NA", "--", "Inf")
  )
  res <- doGSEA(pval_table, gene_list = gene_list,
                tissue_name = tissue_name,
                restrict_ranking = FALSE, 
                top_n = 8)
  
  ggsave(file.path(plot_dir, sprintf("DAS_%s_GSEA.png", t)),
         plot = res$plot, scale = 1, width = 9, height = 8, 
         dpi = 400, bg = "white")
  
  tbls <- make_and_save_tables(
    gsea_res     = res$result,
    tissue_name  = tissue_name,
    padj_thresh  = padj_thresh,
    out_dir      = out_dir
  )
  
  all_tables[[tissue_name]] <- tbls$table %>%
    dplyr::mutate(tissue = tissue_name, .before = 1)
  
}
