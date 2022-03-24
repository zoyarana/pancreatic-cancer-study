# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/Guo2021/")
outdir <- ("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/guo_2021/csvs/")

# libraries
library(readr)
library(dplyr)
library(openxlsx)

# data
patients_survival <- read.xlsx("42003_2021_2557_MOESM7_ESM.xlsx", sheet = "Figure3. Survival")
patients_survival <- mutate(patients_survival, survival_time = Days/30)
colnames(patients_survival)[3] <- "OS"

patients_du12 <- filter(patients_survival, OS == 2 & survival_time < 12)
patients_da24 <- filter(patients_survival, OS == 2 & survival_time >= 24 
                        | OS == 1 & survival_time >= 24)

patients_survival <- patients_survival[patients_survival$Sample_ID %in% patients_da24$Sample_ID 
                                       | patients_survival$Sample_ID %in% patients_du12$Sample_ID,]
patients_survival <- mutate(patients_survival, 
                            status = case_when(patients_survival$OS == 2 
                                               & patients_survival$survival_time < 12 ~ "short_term_survivors", 
                                               patients_survival$OS == 2 & patients_survival$survival_time >= 24 
                                               | patients_survival$OS == 1 & patients_survival$survival_time >= 24 
                                               ~ "long_term_survivors"))
write.csv(patients_survival, paste0(outdir, "patients_survival_st12_lt24.csv"))