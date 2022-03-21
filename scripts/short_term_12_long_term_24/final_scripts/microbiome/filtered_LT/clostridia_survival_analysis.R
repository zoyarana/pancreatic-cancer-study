# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/summarized_counts/")
outdir <- ("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/microbiome/filtered_LT/")

# libraries
library(readr)
library(dplyr)
library("survival")
library("survminer")

# data
PAAD_lcpm <- read.csv("PAAD_lcpm_microbiome_class.csv")
clostridia_lcpm <- PAAD_lcpm[,c("submitter_id", "Clostridia")]

patients <- read.csv("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/patients_survival_st12_lt24.csv", row.names = 1)

patients_lcpm <- patients %>%
  left_join(clostridia_lcpm, by = c("X_PATIENT" = "submitter_id"))

clostridia_median <- median(patients_lcpm$Clostridia)
clostridia_class <- ifelse(patients_lcpm$Clostridia >= clostridia_median, "high", "low")
patients_lcpm <- cbind(patients_lcpm, clostridia_class)

# survival fit
fit <- survfit(Surv(survival_time, OS) ~ clostridia_class, data = patients_lcpm)
summary(fit)

df <- data.frame( time = fit$time, 
                  n.risk = fit$n.risk,
                  n.event = fit$n.event,
                  n.censor = fit$n.censor,
                  surv = fit$surv,
                  upper = fit$upper,
                  lower = fit$lower)

# survival curve
pdf(paste0(outdir,"clostridia_survival_analysis_plot.pdf"))
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"))
dev.off()
