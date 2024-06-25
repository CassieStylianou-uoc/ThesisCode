library(readxl)
CCV <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/RiskCNV_overlapWithCCV/Omara_CCV_original.xlsx",sheet =3) 
CCV <- as.data.frame(CCV)

candGenes <- read.delim("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/UniqueCandidateGenes_0.01.txt")
library(GenomicRanges)
CCV_gr <- makeGRangesFromDataFrame(CCV,keep.extra.columns = T,seqnames.field = "Chr",start.field = "start",end.field = "stop")


risk <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/RiskCNV_overlapWithCCV/RiskAssociated.rds")

risk$chr <- paste("chr",risk$chr,sep="")
risk_sub <- risk[,-c(1)]
colnames(risk_sub) <- c("Chr","Start","Stop","Type")


RA_1_gr <- makeGRangesFromDataFrame(risk_sub,keep.extra.columns = T,seqnames.field = "Chr",start.field = "Start",end.field = "Stop")

mergeByOverlaps(CCV_gr,RA_1_gr)
mergeByOverlaps(RA_1_gr,CCV_gr)

# Calculate distances between CNVs and SNPs
distances <- distanceToNearest(RA_1_gr, CCV_gr)

# Extract the minimum and maximum distances for each SNP
min_distances <- sapply(as.data.frame(distances), function(x) min(x, na.rm = TRUE))
max_distances <- sapply(as.data.frame(distances), function(x) max(x, na.rm = TRUE))

# Convert distances from base pairs to kilobases (if necessary)
min_distances_kb <- min_distances / 1000
max_distances_kb <- max_distances / 1000


# Combine the results into a data frame
summary_df <- data.frame(
  SNP = names(min_distances),
  Min_Distance_KB = min_distances_kb,
  Max_Distance_KB = max_distances_kb
)



############## Testing overlap of all 3+ CNVs to CCVs and seeing if there any overlap #######################

c <-  readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_3_200.rds")
c$chr <- paste("chr",c$chr,sep="")
c$CNV_ID_full <- paste(c$chr,c$startpos,c$endpos,c$cnv_call,sep="_")
c_gr <- makeGRangesFromDataFrame(c,keep.extra.columns = T,seqnames.field = "chr",start.field = "startpos",end.field = "endpos")



CCV_vs_C <- mergeByOverlaps(CCV_gr,c_gr)

summary_df <- as.data.frame(unique(CCV_vs_C$CNV_ID_full))
colnames(summary_df) <- c("CNV_ID")
for(i in 1:nrow(summary_df)){
  summary_df$num.CCVs.overlapped[i] <- length(unique(CCV_vs_C$`Candidate SNP`[CCV_vs_C$CNV_ID_full == summary_df$CNV_ID[i]]))
}

anno <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_1/Gene_Exon_Pli_RegRegions_2probeAnnotation.rds")
anno_ccv <- anno[anno$CNV_ID %in% summary_df$CNV_ID,]
# Sort both data frames by CNV_ID
summary_df_order <- summary_df[order(summary_df$CNV_ID), ]
anno_ccv_order <- anno_ccv[order(anno_ccv$CNV_ID), ]


combined_df <- cbind(summary_df_order, anno_ccv_order)

write.csv(combined_df,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/RiskCNV_overlapWithCCV/ThreeProbe_anno_CCV.csv")
