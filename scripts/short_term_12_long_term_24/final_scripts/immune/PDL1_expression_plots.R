# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/immune/specific_gene_expression_plots/"

# libraries 
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(openxlsx)

# subsetting gene
de <- read.csv("differential_gene_expression_analysis.csv", row.names = 1)
gene <- de %>%
  filter(ENSEMBL == "ENSG00000120217")

# lcpm dataframe
lcpm <- read.csv("lcpm_df.csv", row.names = 1)
lcpm_PDL1 <- lcpm[, c("ENSG00000120217", "status")]

# lcpm vs status boxplot
lcpm_PDL1$status <- as.factor(lcpm_PDL1$status)
xlabs <- paste(levels(lcpm_PDL1$status),"\n(N=",table(lcpm_PDL1[,"status"]),")",sep="")
ggplot(lcpm_PDL1, aes_string(x= "status", y = "ENSG00000120217", color = "status")) + 
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.05)) + stat_compare_means(method = "t.test") +
  theme_bw() +
  scale_x_discrete(labels= xlabs)
ggsave(paste0(outdir, "PDL1_logCPM_boxplot.pdf"))

# lcpm vs IFN gamma response linear regression
gamma_resp <- read.xlsx("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/NIHMS958212-supplement-2 (1).xlsx")
gamma_resp <- gamma_resp[, c(1,13)]

lcpm_PDL1 <- lcpm[, c("ENSG00000120217", "X_PATIENT")]  

PDL1_gamma_resp <- left_join(lcpm_PDL1, gamma_resp, by = c("X_PATIENT" = "TCGA.Participant.Barcode"))

ggscatter(PDL1_gamma_resp, x = "IFN-gamma.Response", y = "ENSG00000120217", 
          add = "reg.line") + 
  stat_cor() 
ggsave(paste0(outdir, "PDL1_gamma_resp_plot.pdf"))
