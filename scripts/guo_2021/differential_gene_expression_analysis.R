# working directories
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/guo_2021/csvs/")
# outdir <- ("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/guo_2021/csvs/")

# libraries
library(readr)
library(limma)
library(edgeR)
library(dplyr)
library(Homo.sapiens)
library(tibble)

# data
patients_survival <- read.csv("patients_survival_st12_lt24.csv", row.names = 1)
pid_patients <- patients_survival$Sample_ID

#  counts
counts <- read_tsv("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/Guo2021/GSE172356_PDA_gene_expression_matrix.txt")
counts <- column_to_rownames(counts, "gene")
pid_counts <- colnames(counts)

pid <- match(pid_patients, pid_counts)
counts <- counts[,pid]
counts <- na.omit(counts)

genes <- select(Homo.sapiens, keys=rownames(counts), columns=c("SYMBOL", "TXCHROM"),
                keytype="SYMBOL")
genes <- genes[!duplicated(genes$SYMBOL),]

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
status<- as.factor(patients_survival[,5])

saveRDS(dge, "dge.RDS")

# lcpm data frame
lcpm_df <- data.frame(lcpm)
lcpm_df <- data.frame(t(lcpm_df))
lcpm_df <- left_join(rownames_to_column(lcpm_df), patients_survival, by = c("rowname" = "Sample_ID"))
colnames(lcpm_df)[1] <- "Sample_ID"
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
  mutate(SYMBOL  = rownames(st.vs.lt)) %>%
  left_join(genes)

# top genes filtering and labels
st.vs.lt <- mutate(st.vs.lt, score = abs(logFC)*(1-adj.P.Val)) %>%
  arrange(st.vs.lt, desc(score))
write.csv(st.vs.lt,"differential_gene_expression_analysis.csv")
