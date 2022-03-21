# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/immune/RNA_signature/"

# libraries
library(readr)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(openxlsx)

# subset patients
patients_survival <- read.csv("patients_survival_st12_lt24.csv", row.names= 1)

# Add leukocyte fraction
infiltrate <- read.xlsx("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/NIHMS958212-supplement-2 (1).xlsx")
infiltrate <- infiltrate[,!duplicated(colnames(infiltrate))]

patients_survival$BARCODE <- gsub("(TCGA\\-[[:alnum:]]{2}\\-[[:alnum:]]{4}).*","\\1",
                                  patients_survival$SAMPLE_ID)

annot <- patients_survival %>%
  left_join(infiltrate, by = c("BARCODE" = "TCGA.Participant.Barcode")) 

# boxplots
colnames(annot) <- gsub("-","",colnames(annot))
annot[,"status"] <- factor(annot[,"status"])

for(col in 16:73) {
  colname <- colnames(annot)[col]
  xlabs <- paste(levels(annot[,"status"]),"\n(N=",table(annot[,"status"]),")",sep="")
  ggplot(annot, aes_string(x= "status", y = colname, color = "status")) + geom_boxplot() +
    geom_point(position = position_jitter(w = 0.05)) + stat_compare_means(method = "t.test") +
    theme_bw() +
    scale_x_discrete(labels=xlabs)
  ggsave(paste0(outdir, colname,"_boxplot.pdf"))
}

