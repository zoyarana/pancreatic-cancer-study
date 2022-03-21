# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/dge/normalisation_plots/"

# libraries
library(readr)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggrepel)
library(RColorBrewer)
library(tibble)
library(pheatmap)

# data
patients_survival <- read.csv("patients_survival_st12_lt24.csv", row.names = 1)

dge <- readRDS("dge.RDS")

# logCPM density
nsamples <- ncol(dge)
col <- brewer.pal(nsamples, "Paired")

cpm <- cpm(dge)
lcpmu <- cpm(dge, log = TRUE)

pdf(paste0(outdir, "filtered_genes.pdf"))
par(mfrow=c(1,2))
plot(density(lcpmu[,5]), col=col[1], lwd=2, ylim=c(0,0.5), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpmu[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
lcpmf <- cpm(dge, log=TRUE)
plot(density(lcpmf[,5]), col=col[1], lwd=2, ylim=c(0,0.2), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
  den <- density(lcpmf[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
dev.off()

# TMM method 
dge2 <- dge

pdf(paste0(outdir,"TMM_boxplots.pdf"))
par(mfrow=c(1,2))
lcpm <- cpm(dge, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Unnormalised data",ylab="Log-cpm")

dge2 <- calcNormFactors(dge2, method = "TMM")
dge2$samples$norm.factors

lcpm <- cpm(dge2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Normalised data",ylab="Log-cpm")
dev.off()

# MDS plot
# MDS based on long-term vs. short-term survivors
status<- as.factor(patients_survival[,11])
lcpm <- cpm(dge, log=TRUE)
pdf(paste0(outdir,"MDS_plot_LtvsSt.pdf"))
par(mfrow=c(1,1))
col.status <- status
levels(col.status) <-c("red","black")
col.status <- as.character(col.status)
plotMDS(lcpm, labels=status, col=col.status )
title(main="Patient Survival Status")
dev.off()

# MDS based on survival time
time<- as.factor(round(patients_survival[,"survival_time"], 0))
pdf(paste0(outdir,"MDS_plot_survivaltime.pdf"))
par(mfrow=c(1,1))
col.status <- status
levels(col.status) <-c("red","black")
col.status <- as.character(col.status)
plotMDS(lcpm, labels=time, col=col.status )
title(main="Patient Survival Status")
dev.off()

# MDS based on patient id
patient_id <- patients_survival$SAMPLE_EXP
pdf(paste0(outdir,"MDS_plot_patientid.pdf"))
par(mfrow=c(1,1))
col.status <- status
levels(col.status) <-c("red","black")
col.status <- as.character(col.status)
plotMDS(lcpm, labels= patient_id, col=col.status )
title(main="Patient Survival Status")
dev.off()

# volcano plot
st.vs.lt <- read.csv("differential_gene_expression_analysis.csv", row.names = 1)

topgenes <- slice_max(st.vs.lt, order_by = score, n = 10)
toplabels <- pull(topgenes, SYMBOL)    
st.vs.lt <- mutate(st.vs.lt, SYMBOL = ifelse(SYMBOL %in% toplabels, SYMBOL, ""))

# significance and colors
down <- unlist(strsplit('Down,Not Sig,Up', split=","))[1]
notsig <- unlist(strsplit('Down,Not Sig,Up', split=","))[2]
up <- unlist(strsplit('Down,Not Sig,Up', split=","))[3]
st.vs.lt <- mutate(st.vs.lt, sig = case_when(
  1-adj.P.Val > 0.95 & logFC > 1 ~ up, 
  1-adj.P.Val > 0.95 & logFC < -1 ~ down,
  TRUE ~ notsig))
pdf(paste0(outdir,"volcano_plot.pdf"))
ggplot(st.vs.lt, aes(x= logFC, y= 1-adj.P.Val))+
  geom_point(aes(color = sig)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line("black"),
        legend.key = element_blank()) +
  scale_colour_manual(values = c("blue","grey","red")) +
  geom_text_repel(data = filter(st.vs.lt, SYMBOL!=""), aes(label=SYMBOL),
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  show.legend = FALSE)
dev.off()


