# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/immune/digital_path/"

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
patients_survival <- read.csv("patients_survival_st12_lt24.csv", row.names = 1)

# add cell types data
cell_types <- read.xlsx("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/TILMap_TableS1_Thorsson_digital_path.xlsx")

patients_survival <- patients_survival %>% 
  left_join(cell_types, by = c("X_PATIENT" = "ParticipantBarcode"))

# survival time v. cell types linear regression
for(col in 17:22) {
  colname <- colnames(patients_survival)[col]
  ggscatter(patients_survival, x = "survival_time", y = colname, 
            add = "reg.line") + 
    stat_cor() 
  ggsave(paste0(outdir, colname,"_linreg_plot.pdf"))
}

# status v.cell type boxplots
patients_survival[,"status"] <- factor(patients_survival[,"status"])

for(col in 17:22) {
  colname <- colnames(patients_survival)[col]
  xlabs <- paste(levels(patients_survival[,"status"]),"\n(N=",table(patients_survival[,"status"]),")",sep="")
  ggplot(patients_survival, aes_string(x= "status", y = colname, color = "status")) + geom_boxplot() +
    geom_point(position = position_jitter(w = 0.05)) + stat_compare_means(method = "t.test") +
    theme_bw() +
    scale_x_discrete(labels=xlabs)
  ggsave(paste0(outdir, colname,"_boxplot.pdf"))
}
