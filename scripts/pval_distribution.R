library(tidyverse)
library(patchwork)
library(fs)

input_dir <- "/Users/mingruili/Library/CloudStorage/OneDrive-Emory/shiny/data/p_value"
csv_files <- dir_ls(input_dir, regexp = "\\.csv$")
data_list <- map(csv_files, ~read_csv(.x) %>% mutate(filename = path_file(.x)))

# 1. Create cleaner individual plots
plot_list <- imap(data_list, ~{
  ggplot(.x, aes(x = p_value_age)) + # Remove as.numeric() if already numeric
    geom_histogram(bins = 15, fill = "#1f77b4", alpha = 0.9) + # Fewer bins
    labs(title = str_remove(path_file(.y), "_pvalue_results\\.csv$")) +
    # theme_void() + # Minimal theme
    theme(
      plot.title = element_text(size = 4, margin = margin(b = 2)),
      plot.margin = margin(1, 1, 1, 1, "mm")
    )
})

# 2. Optimized combination and saving
combined_plot <- wrap_plots(plot_list, ncol = 9, nrow = 6) +
  plot_annotation(title = "raw p-value Distributions") &
  theme(plot.title = element_text(hjust = 0.5, size = 12))

# 3. Save with sane dimensions
ggsave("combined_plots.png", combined_plot,
       path = input_dir,
       width = 16, height = 12, # Reduced from 36x24
       dpi = 400, # Higher DPI instead of huge dimensions
       limitsize = FALSE)
