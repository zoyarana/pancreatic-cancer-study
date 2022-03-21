# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/immune/specific_gene_expression_plots/"

# libraries 
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# subsetting gene
de <- read.csv("differential_gene_expression_analysis.csv", row.names = 1)
gene <- de %>%
  filter(ENSEMBL == "ENSG00000115415")

# lcpm dataframe
lcpm <- read.csv("lcpm_df.csv", row.names = 1)
lcpm_stat1 <- lcpm[, c("ENSG00000115415", "status")]

# lcpm vs status
lcpm_stat1$status <- as.factor(lcpm_stat1$status)
xlabs <- paste(levels(lcpm_stat1$status),"\n(N=",table(lcpm_stat1[,"status"]),")",sep="")
ggplot(lcpm_stat1, aes_string(x= "status", y = "ENSG00000115415", color = "status")) + 
  geom_boxplot() +
  geom_point(position = position_jitter(w = 0.05)) + stat_compare_means(method = "t.test") +
  theme_bw() +
  scale_x_discrete(labels= xlabs)
ggsave(paste0(outdir, "STAT1_logCPM_boxplot.pdf"))

# STAT1 volcano plot
keygene <- c("ENSG00000115415", "ENSG00000120217", "ENSG00000157873", "ENSG00000168961",
             "ENSG00000204287", "ENSG00000114013", "ENSG00000185745", "ENSG00000157601")
genes <- de %>%
  filter(ENSEMBL %in% keygene)
labgene <- pull(genes, SYMBOL)
de <- mutate(de, SYMBOL = ifelse(SYMBOL %in% labgene, SYMBOL, ""))

down <- unlist(strsplit('Down,Not Sig,Up', split=","))[1]
notsig <- unlist(strsplit('Down,Not Sig,Up', split=","))[2]
up <- unlist(strsplit('Down,Not Sig,Up', split=","))[3]
de <- mutate(de, sig = case_when(
  1-adj.P.Val > 0.95 & logFC > 1 ~ up, 
  1-adj.P.Val > 0.95 & logFC < -1 ~ down,
  TRUE ~ notsig))

ggplot(de, aes(x= logFC, y= 1-adj.P.Val))+
  geom_point(aes(color = sig)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line("black"),
        legend.key = element_blank()) +
  scale_colour_manual(values = c("blue","grey","red")) +
  geom_text_repel(data = filter(de, SYMBOL!=""), aes(label=SYMBOL),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  show.legend = FALSE)
ggsave(paste0(outdir, "STAT1_ligands_volcano_plot.pdf"))
