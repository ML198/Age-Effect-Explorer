library(ggplot2)

pval <- read.csv("../data/p_value/cervix_ectocervix_pvalue_results.csv")
pval <- read.csv("../data/p_value/colon_transverse_pvalue_results.csv")

ggplot(pval, aes(x = as.numeric(BH_adjusted))) + 
  geom_histogram(bins = 30, fill = "blue", color = "black") +
  

str(pval$BH_adjusted)
