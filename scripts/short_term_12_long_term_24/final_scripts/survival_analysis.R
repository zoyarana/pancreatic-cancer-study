# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/csvs")

# libraries
library(readr)
library(dplyr)
library("survival")
library("survminer")

# patients
patients_survival <- read.csv("patients_survival_st12.csv", row.names = 1)

# survival fit
fit <- survfit(Surv(survival_time, OS) ~ status, data = patients_survival)
summary(fit)

df <- data.frame( time = fit$time, 
                  n.risk = fit$n.risk,
                  n.event = fit$n.event,
                  n.censor = fit$n.censor,
                  surv = fit$surv,
                  upper = fit$upper,
                  lower = fit$lower)

# survival curve
pdf("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/survival_analysis_plot.pdf")
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"))
dev.off()
