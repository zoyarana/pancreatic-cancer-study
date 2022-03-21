# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/"

# libraries
library(readr)
library(limma)
library(edgeR)
library(dplyr)
library(tibble)

# subset patients
patients <- read.csv("PAAD_annot.csv")
survival <- read_tsv("TCGA-PAAD.survival.tsv")
survival <- mutate(survival, survival_time = OS.time/30)
patients_survival <- inner_join(patients, survival, by = c("SAMPLE_ID" = "sample"), keep = TRUE)

tumor_present <- patients_survival$SAMPLE_EXP[patients_survival$SAMPLE_TYPE == "Tumor"]
patients_survival <- patients_survival[patients_survival$SAMPLE_EXP %in% tumor_present,]

patients_du12 <- filter(patients_survival, OS == 1 & survival_time < 12)
patients_da24 <- filter(patients_survival, OS == 1 & survival_time >= 24 
                        | OS == 0 & survival_time >= 24)

patients_survival <- patients_survival[patients_survival$SAMPLE_EXP %in% patients_da24$SAMPLE_EXP 
                                       | patients_survival$SAMPLE_EXP %in% patients_du12$SAMPLE_EXP,]
patients_survival <- mutate(patients_survival, 
                            status = case_when(patients_survival$OS == 1 
                                               & patients_survival$survival_time < 12 ~ "short_term_survivors", 
                                               patients_survival$OS == 1 & patients_survival$survival_time >= 24 
                                               | patients_survival$OS == 0 & patients_survival$survival_time >= 24 
                                               ~ "long_term_survivors"))
write.csv(patients_survival, paste0(outdir, "patients_survival_st12_lt24.csv"))


