library("hpar")
library("ExperimentHub")
library(tidyverse)
library(ggplot2)
library(ggforce)
library(data.table)
library(tidyverse)
library(cBioPortalData)
library(cbioportalR)
library(dplyr)
library(readxl)
library(purrr)
library(AnVIL)

'%!in%' <- function(x,y)!('%in%'(x,y))
############### HUMAN PROTEIN ATLAS - ENDOMETRIAL TISSUE 

hpa_data <- allHparData()
data(rnaGeneTissue)
head(rnaGeneTissue)

rna <- rnaGeneTissue
rna_endometrial <- rna[rna$Tissue == "endometrium",]
write.csv(rna_endometrial,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/Full_HPA_endometrium.csv")
candidatGenes <- read.table("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/UniqueCandidateGenes_0.01.txt",header = T)

HPA_wanted <- rna[rna$Gene.name %in% candidatGenes$Genes & rna$Tissue== "endometrium",]
length(unique(HPA_wanted$Gene.name))
length(unique(HPA_wanted$Gene.name)) #141

# Define cutoffs for different categories
cutoffs <- c(0.5, 10, 1000)

HPA_wanted_2 <- HPA_wanted

HPA_wanted_2$Expression_category <- cut(HPA_wanted_2$TPM, breaks = c(0, cutoffs, Inf), labels = c("Below Cutoff", "Low Expression", "Medium Expression", "High Expression"), right = FALSE)



ExpressedInNormal <- HPA_wanted_2[HPA_wanted_2$Expression_category %in% c("Low Expression","Medium Expression"),]
### ALDOA has two TMSB15B has two - picked these: ENSG00000158427 + ENSG00000149925 (removed ENSG00000285043,ENSG00000269226)

ExpressedInNormal<- ExpressedInNormal[-c(which(ExpressedInNormal$Gene %in% c("ENSG00000285043","ENSG00000269226"))),]

candidatGenes_expressed <- ExpressedInNormal$Gene.name



############ DOSAGE SENSITIVITY
set_cbioportal_db(db = "public")

u_panCancer <- cBioDataPack("ucec_tcga_pan_can_atlas_2018",ask=FALSE)
u_panCancer ## see whats available 
CNV <- assay(u_panCancer[[2]])
mRNA <- assay(u_panCancer[[7]])

write.csv(CNV,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/Full_CNA_endoTumour.csv")
write.csv(CNV,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/Full_mRNA_endoTumour.csv")

CNV_want <- CNV[which(rownames(CNV)%in% candidatGenes_expressed),] #107 genes / 523 samples
mRNA_want <- mRNA[which(rownames(mRNA) %in% candidatGenes_expressed),] #107 genes / 527 samples

common_cols <- intersect(colnames(CNV_want),colnames(mRNA_want))
common_rows <- intersect(rownames(CNV_want),rownames(mRNA_want))
CN_common <- CNV_want[common_rows,common_cols] #107 / 521
Exp_common <- mRNA_want[common_rows,common_cols] #107 / 521


setdiff(candidatGenes_expressed,rownames(CN_common))


########### Some genes have no data present 
t <- as.data.frame(Exp_common)

for(i in 1:nrow(t)){
  t$empty.count[i] <- sum(is.na(t[i,]))
}
table(t$empty.count)

LackingExpData <- rownames(t)[t$empty.count>=1]
## No expression data for these genes: [1] "CHRM2" "CRHR1"  

Exp_edit <- Exp_common[-c(which(rownames(Exp_common) %in% LackingExpData)),] # 105 521 
CN_edit <- CN_common[c(rownames(CN_common) %in% rownames(Exp_edit)),] #105 521




#########################
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

## set p-value threshold
pthreshold = 0.0001


# Replace all occurrences of -2 with -1 across all columns
CN_edit[CN_edit == -2] <- -1
CN_edit[CN_edit == 2] <- 1


# Function to summarize copy numbers for each gene
summary_by_gene <- function(df) {
  summary <- apply(df[,], 2, function(x) table(x))
  return(summary)
}


# Apply the function to get the summary
gene_summary <- summary_by_gene(t(CN_edit))

##### Given there are very few (deep deletion - indicative of homozygous??) - decided to collapse -2 and -1 into "deletion" category 



#############################################################################
## data wrangling for plotting
#############################################################################
## Generate data frame for plotting
RNA <- as.data.frame(t(Exp_edit))
CNA <- as.data.frame(t(CN_edit))
#CNA <- as.data.frame(CNA_reformat)
#CNA <- as.data.frame(t(CN_edit))
p_df <- cbind(gather(RNA),gather(CNA))
colnames(p_df)<- c("gene", "RNA", "gene2", "CN")
p_df$CN <- factor(p_df$CN)
p_df$gene <- factor(p_df$gene)

library(dplyr)


# Assuming your dataframe is named df
df_mean <- p_df %>%
  group_by(gene, CN) %>%
  summarize(mean_RNA = mean(RNA, na.rm = TRUE),
            frequency = n()) %>%
  pivot_wider(names_from = CN, values_from = c(mean_RNA, frequency), 
              names_sep = "_", names_prefix = "CN_")


write.csv(df_mean,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/DosageSensitivity/ZscoreMrna_CN_TCGA_endo.csv")

## linear models for gene dosage
models <- lapply(split(p_df, p_df$gene), function(i) lm(RNA~CN, data=i)) #105


# Initialize a vector to store the indices of models with negative coefficients
models_with_negative_coefficients <- numeric()


# Iterate over each model in the list
for (i in seq_along(models)) {
  # Extract coefficients for CN
  coefficients <- coef(models[[i]])
  
  # Check if any coefficient for CN is negative (excluding intercept)
  if (any(coefficients[-1] < 0)) {
    models_with_negative_coefficients <- c(models_with_negative_coefficients, i)
  }
}

# Print the indices of models with negative coefficients
if (length(models_with_negative_coefficients) > 0) {
  cat("Models with negative coefficients for CN:", models_with_negative_coefficients, "\n")
} else {
  print("No models have negative coefficients for CN.")
} # Models with negative coefficients for CN: 6 9 18 23 32 41 57 59 61 63 76 79 86 94 100 101 102 


reshape_all.pvalues <- lapply(models, function(i) lmp(i))
reshape_all.pvalues <- melt(data.frame(reshape_all.pvalues))

models_positive <- models[-models_with_negative_coefficients] # 88 remain



# Extract p-values and make a geom_text df
library(reshape)
p.values <- lapply(models_positive, function(i) lmp(i))
p.values <- melt(data.frame(p.values))
colnames(p.values) <- c('gene', 'pvalue')
p.values$pvalue2 <- format(p.values$pvalue, digits = 2)

p_df_reformat <- p_df
p_df_reformat <- p_df_reformat[p_df_reformat$gene %in% p.values$gene,]
p_df_reformat_df <- as.data.frame(p_df_reformat)



# Find the maximum RNA value for each gene in p_df
max_RNA <- aggregate(RNA ~ gene, data = p_df_reformat, FUN = max)

# Merge the maximum RNA values into the p.values dataframe based on gene names
p.values <- merge(p.values, max_RNA, by = "gene", all.x = TRUE)
# Rename the merged column to y_max
colnames(p.values)[4] <- "ymax"
signif.gene <- p.values$gene[p.values$pvalue < pthreshold] # 59
p.values$pvalue2 <- format(p.values$pvalue, digits = 2)
length(signif.gene) #59



DS.plots <- ggplot(p_df_reformat[p_df_reformat$gene %in% signif.gene,],aes(x = CN, y = RNA)) + geom_boxplot(fill = "#B0C4De",outlier.size = 0.75, size=.1) +  theme_bw()+
  geom_smooth(method = "lm",se=FALSE, col = "black", aes(group=1), size=0.75)+facet_wrap_paginate(~gene, scales='free', ncol=4, nrow=5)+
  geom_text(data = p.values[p.values$gene %in% signif.gene,], aes(x=1.2,y=ymax-0.5,label=pvalue2),size = 3)



DS.plots <- ggplot(p_df_reformat[p_df_reformat$gene %in% signif.gene,],aes(x = CN, y = RNA)) + geom_boxplot(fill = "#B0C4De",outlier.size = 0.75, size=.1) +  theme_bw()+
  geom_smooth(method = "lm",se=FALSE, col = "black", aes(group=1), size=0.75)+facet_wrap_paginate(~gene, scales='free', ncol=4, nrow=5)+
  geom_text(data = p.values[p.values$gene %in% signif.gene,], aes(x=1.2,y=ymax-0.5,label=paste0("p = ",pvalue2)),size = 3)+  labs(x = "Gene copy number (GISTIC)", y = "Candidate gene RNA expression (z-scores)")


pdf("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/DosageSensitivity/DosageSensitivityOf59genes_axis.pdf",paper='a4',width = 5 * 8.27, height = 8 * 11.69)
for(i in 1:n_pages(DS.plots)){
  print(DS.plots+facet_wrap_paginate(~gene, scales='free', ncol=4, nrow=5, page=i))
}
dev.off()







################# marrying up expression and dosage sensitivity data for supplementary info
CG_normal_expres <-read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/CandidateGene_141_ExpressionNormalEndo.csv")

for(i in 1:nrow(CG_normal_expres)){
  if(CG_normal_expres$Gene.name[i] %in% p.values$gene){
    CG_normal_expres$DS.p[i] <- p.values$pvalue2[p.values$gene == CG_normal_expres$Gene.name[i]]
  }else{
    CG_normal_expres$DS.p[i] <- c("N/A")
  }
}

###############################################















