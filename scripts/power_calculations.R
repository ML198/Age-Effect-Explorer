library(tidyverse)
library(pwr)

# Load sample sizes data file
ss <- tissue_sizes_df

# Named vector: names are tissues, values are sample sizes
tissue_ss <- setNames(ss$Sample_size, ss$Tissue)

# Significance parameters
alpha <- 0.05           # significance level
k_preds <- 1            # number of predictors of interest (1 because of 'age' predictor)

# Cohen's f^2 effect sizes
f2_values <- c(small = 0.02, medium = 0.15, large = 0.35)

# Compute power for each tissue and effect size
results_list <- lapply(seq_along(tissue_ss), function(i) {
  
  tissue_name = names(tissue_ss)[i]
  n <- tissue_ss[i]
  
  # denominator df v = N - k_preds - 1
  v <- n - k_preds - 1
  
  # handle too-small sample sizes
  if (v <= 0) {
    # Not enough observations to estimate; return NA row
    return(data.frame(
      tissue = tissue_name,
      sample_size = n,
      k_preds = k_preds,
      v = v,
      f2_label = names(f2_values),
      f2 = as.numeric(f2_values),
      power = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  # compute power for each effect size
  powers <- sapply(f2_values, function(f2) {
    pwr_res <- pwr.f2.test(u = k_preds, v = v, f2 = f2, sig.level = alpha)
    # pwr.f2.test returns $power
    pwr_res$power
  })
  
  df <- data.frame(
    tissue = tissue_name,
    sample_size = n,
    k_preds = k_preds,
    v = v,
    effect_size = names(f2_values),
    f2 = as.numeric(f2_values),
    power = as.numeric(powers),
    stringsAsFactors = FALSE
  )
  return(df)
})

results_df <- Reduce(rbind, results_list) %>%
  arrange(tissue, f2)

write.csv(results_df, "data/summary/power_calculations.csv", col.names = FALSE)
res <- list(sample_sizes = ss, power_df = results_df)
saveRDS(res, "data/summary/power_calculations.rds")



library(tidyverse)

power_calculations[["power_df"]] %>%
  ggplot(aes(x = sample_size, y = power, color = effect_size)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Statistical Power by Sample Size across GTEx Tissues",
    subtitle = "Effect sizes: small (f²=0.02), medium (f²=0.15), large (f²=0.35)",
    x = "Sample size per tissue",
    y = "Power"
  )

power_calculations[["power_df"]] %>%
  ggplot(aes(x = effect_size, y = reorder(tissue, sample_size), fill = power)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Heatmap of Statistical Power across Tissues and Effect Sizes",
    x = "Effect Size",
    y = "Tissue"
  )
