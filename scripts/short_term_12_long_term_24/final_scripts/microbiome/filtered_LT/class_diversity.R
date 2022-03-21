# working directory
setwd("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/data/summarized_counts/")
outdir <- ("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/microbiome/filtered_LT/")

# libraries
library(readr)
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(tibble)
library(phyloseq)
library(car)
library(ggpubr)

# data
class_abundance <- read.csv("PAAD_lcpm_abundance_microbiome_class.csv")

patients <- read.csv("/bioinf_core/UROP/Zoya_Rana/01_pancreatic_cancer/results/short_term_12_long_term_24/final_results/csvs/patients_survival_st12_lt24.csv", row.names = 1)
patients <- patients %>%
  filter(SAMPLE_EXP == "TCGA.3A.A9IJ.01A.11R.A39D.07"| SAMPLE_EXP == "TCGA.3A.A9IS.01A.21R.A39D.07"| SAMPLE_EXP == "TCGA.3A.A9IO.01A.11R.A38C.07"| 
           SAMPLE_EXP == "TCGA.3A.A9IR.01A.11R.A38C.07"| SAMPLE_EXP == "TCGA.3A.A9IL.01A.11R.A38C.07"| SAMPLE_EXP == "TCGA.2L.AAQM.01A.11R.A39D.07"|
           SAMPLE_EXP == "TCGA.3A.A9IV.01A.11R.A41B.07" | status == "short_term_survivors" )

class_abundance <-class_abundance %>%
  inner_join(patients, by = c("submitter_id" = "X_PATIENT"))

# shannon 
class_abundance$alpha_shan <- diversity(class_abundance[,2:57], 
                                        MARGIN = 1, index = "shannon")
hist(class_abundance$alpha_shan)
qqPlot(lm(alpha_shan ~ status, data = class_abundance), simulate = TRUE, main = "QQ Plot", labels = FALSE)

fit_shannon<-aov(alpha_shan ~ status, data = class_abundance)
TukeyHSD(fit_shannon)
class_abundance$status <- as.factor(class_abundance$status)

pdf(paste0(outdir, "class_shannon_plot.pdf"))
ggboxplot(class_abundance, x="status",y="alpha_shan",color="status",outlier.shape = NA)+
  stat_compare_means(method = "t.test")+
  scale_color_manual(values=c("deepskyblue","gray20"))+
  geom_jitter(aes(color=status),size=0.8)+
  ggtitle("Shannon")+ylab("")+xlab("")+guides(color=F)+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank(),
        panel.background = element_rect(fill = "grey90",colour = NA), panel.border = element_blank(),
        panel.grid.major = element_line(colour = "white"),axis.line = element_blank(),
        legend.position = "right")
dev.off()


# distance and PCoA
distance <- vegdist(class_abundance[,2:57], method = "bray")
distance <- as.matrix(distance)

pcoa <- cmdscale(distance, k = 2, eig = TRUE)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_status <- data.frame({pcoa$point})[1:2]
sample_status$names <- class_abundance$submitter_id
names(sample_status)[1:2] <- c('PCoA1', 'PCoA2')

sample_status <- left_join(sample_status, class_abundance, by= c("names" ="submitter_id"))
colnames(sample_status)[69]<-"group"
sample_status$group <- factor(sample_status$group, levels = c("short_term_survivors", "long_term_survivors"))

pdf(paste0(outdir,"class_PCoA.pdf"))
ggplot(sample_status,aes(PCoA1,PCoA2, color = group))+
  geom_point(aes(color = group), size = 1.5, alpha = 0.7)+
  stat_ellipse(aes( PCoA1,PCoA2,fill=group),geom="polygon",level=0.9,alpha=0.1)+
  scale_fill_manual(values=c("gray20", "deepskyblue"))+
  scale_color_manual(values = c("gray20", "deepskyblue"))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme(legend.title = element_blank())
dev.off()

