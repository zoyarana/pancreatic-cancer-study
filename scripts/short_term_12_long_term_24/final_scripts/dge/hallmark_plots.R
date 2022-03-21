# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/csvs")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/dge/hallmark/"

# libraries
library(readr)
library(gprofiler2)
library(dplyr)
library(ggplot2)
source("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/scripts/short_term_12_long_term_24/pathways/hallmark/hallmark_grouping.R")

# subsetting significant genes
de <- read.csv("differential_gene_expression_analysis.csv", row.names = 1)
genes <- de %>% 
  filter(adj.P.Val < 0.05)

genes_up <- genes %>% filter(logFC >= 0)
genes_down <- genes %>% filter (logFC <0)

# Hallmark – save as csv
# overexpressed genes
query_up <- list(genes_up$SYMBOL)
hallmark_up <- "gp__5vEZ_N4AY_ye0"

gostres_up <- gost(query = query_up,
                   organism = hallmark_up, ordered_query = FALSE,
                   multi_query = FALSE, significant = FALSE, correction_method = "fdr")

results_up <- data.frame(gostres_up$result)
results_up <- group_into_hallmark_categories(results_up)
write.csv(results_up[,-14], paste0(outdir, "hallmark_results_up.csv"))

# underexpressed genes
query_down <- list(genes_down$SYMBOL)
hallmark_down <- "gp__5vEZ_N4AY_ye0"

gostres_down <- gost(query = query_down,
                     organism = hallmark_down, ordered_query = FALSE,
                     multi_query = FALSE, significant = FALSE, correction_method = "fdr")

results_down <- data.frame(gostres_down$result)
results_down <- group_into_hallmark_categories(results_down)
write.csv(results_down[,-14], paste0(outdir, "hallmark_results_down.csv"))


# bar plot
# overexpressed genes
results_up <- results_up %>%
  arrange(ProcessCategory, p_value)
results_up$term_id <- factor(results_up$term_id, levels = results_up$term_id)

pdf(paste0(outdir,"hallmark_up_barplot.pdf"))
ggplot(results_up, aes(x = -log10(p_value), y = term_id, fill = ProcessCategory,
                       group = ProcessCategory)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() +
  geom_vline(aes(xintercept = -log10(0.05), color = "red"))
dev.off()

# underexpressed genes
results_down <- results_down%>%
  arrange(ProcessCategory, p_value)
results_down$term_id <- factor(results_down$term_id, levels = results_down$term_id)

pdf(paste0(outdir,"hallmark_down_barplot.pdf"))
ggplot(results_down, aes(x = -log10(p_value), y = term_id, fill = ProcessCategory)) +
  geom_bar(stat = "identity", color = "black") + 
  theme_minimal() +
  geom_vline(aes(xintercept = -log10(0.05), color = "red"))
dev.off()

