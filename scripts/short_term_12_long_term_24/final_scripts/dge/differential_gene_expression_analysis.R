# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/")
outdir <- "/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/dge/normalisation_plots/"

# libraries
library(readr)
library(limma)
library(edgeR)
library(dplyr)
library(Homo.sapiens)
library(tibble)

# data
patients_survival <- read.csv("patients_survival_st12_lt24.csv", row.names = 1)
pid_patients <- patients_survival$SAMPLE_EXP

#  counts
counts <- read.csv("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/RawCounts.csv", row.names = 1)
pid_counts <- colnames(counts)

pid <- match(pid_patients, pid_counts)
counts <- counts[,pid]

genes <- select(Homo.sapiens, keys=rownames(counts), columns=c("SYMBOL", "TXCHROM"),
                keytype="ENSEMBL")
genes <- genes[!duplicated(genes$ENSEMBL),]

# counts DGE List
dge <- DGEList(counts = counts)

# patient DGE List
dge$patients_survival <- patients_survival

# genes DGE List
dge$genes <- genes

# counts to cpm
cpm <- cpm(dge)

# removing not expressed genes 
# min library size
L <- min(colSums(counts)) 
# 20% population
P <- round(ncol(counts)*0.50) 
keep <- rowSums(cpm(dge) > 5/L*1e6) > P
dge <- DGEList(dge[keep,,keep.lib.sizes=F])
lcpm <- cpm(dge, log = TRUE)
status<- as.factor(patients_survival[,11])

saveRDS(dge, "dge.RDS")

# lcpm data frame
lcpm_df <- data.frame(lcpm)
lcpm_df <- data.frame(t(lcpm_df))
lcpm_df <- left_join(rownames_to_column(lcpm_df), patients_survival, by = c("rowname" = "SAMPLE_EXP"))
write.csv(lcpm_df, "lcpm_df.csv")

# differential gene expression analysis
design <- model.matrix(~0+status)
colnames(design) <- gsub("status", "", colnames(design))
contr.matrix <- makeContrasts(Short_termvsLong_term = 
                                short_term_survivors-long_term_survivors,
                              levels = colnames(design))

vfit <- lmFit(lcpm, design)
vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit, trend = T)
dt <- topTable(efit, p.value = 1)
summary(dt$adj.P.Val < 0.05)

# Volcano plot
st.vs.lt <- topTreat(efit, n = Inf)

# Add genes
st.vs.lt <- st.vs.lt %>%
  mutate(ENSEMBL  = rownames(st.vs.lt)) %>%
  left_join(genes)

# top genes filtering and labels
st.vs.lt <- mutate(st.vs.lt, score = abs(logFC)*(1-adj.P.Val)) %>%
  arrange(st.vs.lt, desc(score))
write.csv(st.vs.lt,"differential_gene_expression_analysis.csv")

# heatmap
dge$genes <- genes
v <- voom(dge, design, plot = TRUE)

st.vs.lt.tg <- st.vs.lt$SYMBOL[1:100]
t <- which (v$genes$SYMBOL %in% st.vs.lt.tg) [1:100]


pdf(paste0(outdir,"heatmap_patients_id.pdf"), height= 7)
heatmap.2(v$E[t,], scale = "row", labRow = v$genes$SYMBOL[t], 
          labCol = patients_survival$SAMPLE_EXP, col = colorpanel(100, "blue", "grey", "red"),
          trace = "none", density.info = "none", margin=c(8,6), lhei = c(2,10),
          dendrogram="column", cexCol=0.5)
dev.off()

pdf(paste0(outdir,"heatmap_status.pdf"), height= 7)
heatmap.2(v$E[t,], scale = "row", labRow = v$genes$SYMBOL[t], 
          labCol = status, col = colorpanel(100, "blue", "grey", "red"),
          trace = "none", density.info = "none", margin=c(8,6), lhei = c(2,10), 
          dendrogram="column")
dev.off()