# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/summarized_counts/")
outdir <- ("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/microbiome/filtered_LT/")

# libraries
library(readr)
library(dplyr)
library(lefser)
library(SummarizedExperiment)

# data 
PAAD_lcpm <- read.csv("PAAD_lcpm_microbiome_class.csv")

patients <- read.csv("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/patients_survival_st12_lt24.csv", row.names = 1)

patients <- patients %>%
  filter(SAMPLE_EXP == "TCGA.3A.A9IJ.01A.11R.A39D.07"| SAMPLE_EXP == "TCGA.3A.A9IS.01A.21R.A39D.07"| SAMPLE_EXP == "TCGA.3A.A9IO.01A.11R.A38C.07"| 
           SAMPLE_EXP == "TCGA.3A.A9IR.01A.11R.A38C.07"| SAMPLE_EXP == "TCGA.3A.A9IL.01A.11R.A38C.07"| SAMPLE_EXP == "TCGA.2L.AAQM.01A.11R.A39D.07"|
           SAMPLE_EXP == "TCGA.3A.A9IV.01A.11R.A41B.07" | status == "short_term_survivors" )
patients <- patients[,c(8,11)]

patients_lcpm <- patients %>%
  left_join(PAAD_lcpm, by = c("X_PATIENT" = "submitter_id"))
patients_lcpm <- data.frame(patients_lcpm)
patients_lcpm$status <- as.factor(patients_lcpm$status)

tmp <-  patients_lcpm[,c(1,2)]
rownames(tmp) <- tmp$X_PATIENT

patients_lcpm_se <- SummarizedExperiment(assays = list(exprs = t(patients_lcpm[,-c(1,2)])),
                                         colData = tmp)
rownames(patients_lcpm_se) <- colnames(patients_lcpm[,-c(1,2)])
colnames(patients_lcpm_se) <- patients_lcpm[,1]
levels(colData(patients_lcpm_se)$status)

# contingency table
table(patients_lcpm$status)

# lefser analysis
res <- lefser(patients_lcpm_se, groupCol = "status", kruskal.threshold = 0.05,
              wilcox.threshold = 0.05,  lda.threshold = 0)

# lefser plot
pdf(paste0(outdir, "class_lda_analysis_0_LT_1_ST_0_05_p.pdf"))
lefserPlot(res)
dev.off()

