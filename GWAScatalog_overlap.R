


library( readxl)
library(GenomicRanges)
dirLift <- "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Original_fromGWAScat/"

EC <- read_excel(paste0(dirLift,"ReadyForLift_GWAScat.xlsx"),sheet =1)
EC <- EC[,c(6:9)]
EC_gr <- makeGRangesFromDataFrame(EC,seqnames.field = "hg19.chr",start.field = "hg19.pos1",end.field = "hg19.pos2",keep.extra.columns = T)

OB <- read_excel(paste0(dirLift,"ReadyForLift_GWAScat.xlsx"),sheet =2)
OB <- OB[,c(6:9)]
OB_gr <- makeGRangesFromDataFrame(OB,seqnames.field = "hg19.chr",start.field = "hg19.pos1",end.field = "hg19.pos2",keep.extra.columns = T)


TD <- read_excel(paste0(dirLift,"ReadyForLift_GWAScat.xlsx"),sheet =3)
TD <- TD[,c(6:9)]
TD$hg19.chr <- gsub("chrX","chr23",TD$hg19.chr)
TD_gr <- makeGRangesFromDataFrame(TD,seqnames.field = "hg19.chr",start.field = "hg19.pos1",end.field = "hg19.pos2",keep.extra.columns = T)


additional_EC <- read_excel(paste0(dirLift,"ReadyForLift_GWAScat.xlsx"),sheet =4)
additional_EC <- additional_EC[,c(6:9)]
additional_EC_gr <- makeGRangesFromDataFrame(additional_EC,seqnames.field = "hg19.chr",start.field = "hg19.pos1",end.field = "hg19.pos2",keep.extra.columns = T)

########## Making UCSC formats for lifted over GWAS cat SNPs that were assessed 

GWAScat_UCSC <- function(gwasCatOb){
  toreturn <- gwasCatOb[,c(1:4)]
  toreturn$Zero <- c("0")
  toreturn$plus <- c("+")
  toreturn <- cbind(toreturn,gwasCatOb[,c(2:3)])
  toreturn$type <- c("128,0,128")
  return(toreturn)
}



EC_ucsc <- GWAScat_UCSC(EC)
OB_ucsc <- GWAScat_UCSC(OB)
TD_ucsc <- GWAScat_UCSC(TD)

EC_ucsc <- EC_ucsc[!duplicated(EC_ucsc$hg19.rsID),]
OB_ucsc <- OB_ucsc[!duplicated(OB_ucsc$hg19.rsID),]
TD_ucsc <- TD_ucsc[!duplicated(TD_ucsc$hg19.rsID),]


nrow(EC_ucsc)
nrow(OB_ucsc)
nrow(TD_ucsc)




#library(xlsx)
#write.xlsx(EC_ucsc, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/LiftedOver_UCSCFormat.xlsx", sheetName="EC", row.names=FALSE)
#write.xlsx(OB_ucsc, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/LiftedOver_UCSCFormat.xlsx", sheetName="OB", append=TRUE, row.names=FALSE)
#write.xlsx(TD_ucsc, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/LiftedOver_UCSCFormat.xlsx", sheetName="TD", append=TRUE, row.names=FALSE)





#write.table(RA_dels_UCSC_final,paste0(outputDir,"UCSCFormat_RiskAssociatedDels.txt"),quote = FALSE,row.names = F,sep = '\t',col.names = F)

x <-readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_geneMergeObjects/AllriskAssociated_riskassociated_0.01_3probeGWAS_grFormat.rds")
x$Chr <- paste("chr",x$Chr,sep="")
gcnv <- makeGRangesFromDataFrame(x,seqnames.field = "Chr",start.field = "Start",end.field = "Stop",keep.extra.columns = T)


mergeGWAScat <- function(riskCNV_gr,GWAScat_gr){
  mergeOb <- mergeByOverlaps(riskCNV_gr,GWAScat_gr)
  a <- list()
  for(i in 1:nrow(mergeOb)){
    a[[i]] <- as.data.frame(mergeOb[i,])
  }
  merge_df <- do.call(rbind,a)
  return(merge_df)
}


Merge_EC <- mergeGWAScat(gcnv,EC_gr)
Merge_TD <- mergeGWAScat(gcnv,TD_gr)


unique(Merge_EC$hg19.rsID) # two in total
length(unique(Merge_TD$hg19.rsID)) #33 
length(unique(Merge_TD$riskCNV_gr.cnvID_wType))


Merge_EC_additional <- mergeGWAScat(gcnv,additional_EC_gr)
unique(Merge_EC_additional$hg19.rsID)

############ How many of the 141 candidate genes had a least one ECAC CNV overlapping a GWAS-catalog risk SNP
length(unique(Merge_EC$cnvID_wType))
length(unique(Merge_TD$cnvID_wType))

length(unique(Merge_EC$Gene.Name)) #HNF1B
length(unique(Merge_TD$Gene.Name))


length(unique(Merge_TD$GWAScat_gr.hg19.rsID[Merge_TD$riskCNV_gr.cnvID_wType == "16_29595483_30156963_del"])) 
length(unique(Merge_TD$Gene.Name[Merge_TD$riskCNV_gr.cnvID_wType == "16_29595483_30156963_del"])) 



########### Summary dataframe 
EC_original <- read_excel("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Original_fromGWAScat/Compiled_2024_GWAScatOriginalData.xlsx",sheet = 3)

sum_EC <- as.data.frame(unique(Merge_EC$hg19.rsID))
colnames(sum_EC) <- c("RiskSNP")

EC_CNVList <- list()
for(i in 1:nrow(sum_EC)){
  sum_EC$Num.Risk.CNVs.Overlapping[i] <- length(unique(Merge_EC$cnvID_wType[Merge_EC$hg19.rsID == sum_EC$RiskSNP[i]]))
  sum_EC$Num.Unique.Cases.WithOverlappingCNV[i] <- length(unique(Merge_EC$sampleid[Merge_EC$hg19.rsID == sum_EC$RiskSNP[i] & Merge_EC$cohort==1]))
  sum_EC$Num.Unique.Controls.WithOverlappingCNV[i] <- length(unique(Merge_EC$sampleid[Merge_EC$hg19.rsID == sum_EC$RiskSNP[i] & Merge_EC$cohort==0]))
  sum_EC$SNP.pos[i] <- EC_original$CHR_POS[EC_original$SNPS == sum_EC$RiskSNP[i]]
  sum_EC$Mapped.Gene[i] <- EC_original$MAPPED_GENE[EC_original$SNPS == sum_EC$RiskSNP[i]]
  sum_EC$Region[i] <- EC_original$REGION[EC_original$SNPS == sum_EC$RiskSNP[i]]
  sum_EC$Context[i] <- EC_original$CONTEXT[EC_original$SNPS == sum_EC$RiskSNP[i]]
  OverlappingCNVs <- unique(Merge_EC$cnvID_wType[Merge_EC$hg19.rsID == sum_EC$RiskSNP[i]])
  EC_CNVList[[i]] <- OverlappingCNVs
  sum_EC$OverlappingCNVs[i] <- paste(OverlappingCNVs,collapse = " ; ")
}

EC_vector <- unlist(EC_CNVList)
length(unique(EC_vector)) #9 

T2_original <- read_excel("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Original_fromGWAScat/Compiled_2024_GWAScatOriginalData.xlsx",sheet = 1)
sum_TD <- as.data.frame(unique(Merge_TD$hg19.rsID))
colnames(sum_TD) <- c("RiskSNP")

TD_CNVList <- list()
for(i in 1:nrow(sum_TD)){
  sum_TD$Num.Risk.CNVs.Overlapping[i] <- length(unique(Merge_TD$cnvID_wType[Merge_TD$hg19.rsID == sum_TD$RiskSNP[i]]))
  sum_TD$Num.Unique.Cases.WithOverlappingCNV[i] <- length(unique(Merge_TD$sampleid[Merge_TD$hg19.rsID == sum_TD$RiskSNP[i] & Merge_TD$cohort==1]))
  sum_TD$Num.Unique.Controls.WithOverlappingCNV[i] <- length(unique(Merge_TD$sampleid[Merge_TD$hg19.rsID == sum_TD$RiskSNP[i] & Merge_TD$cohort==0]))
  sum_TD$SNP.pos[i] <- T2_original$CHR_POS[T2_original$SNPS == sum_TD$RiskSNP[i]]
  sum_TD$Mapped.Gene[i] <- T2_original$MAPPED_GENE[T2_original$SNPS == sum_TD$RiskSNP[i]]
  sum_TD$Region[i] <- T2_original$REGION[T2_original$SNPS == sum_TD$RiskSNP[i]]
  sum_TD$Context[i] <- T2_original$CONTEXT[T2_original$SNPS == sum_TD$RiskSNP[i]]
  OverlappingCNVs <- unique(Merge_TD$cnvID_wType[Merge_TD$hg19.rsID == sum_TD$RiskSNP[i]])
  TD_CNVList[[i]] <- OverlappingCNVs
  sum_TD$OverlappingCNVs[i] <- paste(OverlappingCNVs,collapse = " ; ")
}

TD_vector <- unlist(TD_CNVList)
length(unique(TD_vector)) #78

write.csv(sum_EC,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ECgwasCat_mergedWithRiskCNVs.csv")
write.csv(sum_TD,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/TDgwasCat_mergedWithRiskCNVs.csv")



############## Reading in merge GWAS objects 

Delmerge <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_geneMergeObjects/DelsOnly_3probeGWAS_geneMerge.rds")
Dupmerge <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_geneMergeObjects/DupsOnly_3probeGWAS_geneMerge.rds")
Lofmerge <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_geneMergeObjects/LOF_3probeGWAS_geneMerge.rds")


############ Which candidate genes are involved when we use the GWAS-cat overlapping CNVs

unique(Delmerge$Gene.Name[Delmerge$cnvID_wType %in% EC_vector])
unique(Delmerge$Gene.Name[Delmerge$cnvID_wType %in% TD_vector])



ThreeProbeCNVs <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
ThreeProbeCNVs <- ThreeProbeCNVs[ThreeProbeCNVs$probes>=3,]
ThreeProbeCNVs$cnv_id <- paste(ThreeProbeCNVs$chr,ThreeProbeCNVs$startpos,ThreeProbeCNVs$endpos,ThreeProbeCNVs$cnv_call,sep="_")

#### How many cases and controls carry a CNV that is predicted to overlap a risk SNP

## ENDOMETRIAL CANCER
length(unique(ThreeProbeCNVs$sampleid[ThreeProbeCNVs$cohort==1 & ThreeProbeCNVs$cnv_id %in% EC_vector])) # seven cases
length(unique(ThreeProbeCNVs$sampleid[ThreeProbeCNVs$cohort==0 & ThreeProbeCNVs$cnv_id %in% EC_vector])) # three controls 


## DIABETES
length(unique(ThreeProbeCNVs$sampleid[ThreeProbeCNVs$cohort==1 & ThreeProbeCNVs$cnv_id %in% TD_vector])) # 58 cases
length(unique(ThreeProbeCNVs$sampleid[ThreeProbeCNVs$cohort==0 & ThreeProbeCNVs$cnv_id %in% TD_vector])) # 56 controls 


library(xlsx)
write.xlsx(sum_EC, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Final_RiskSNPOverlap_Corrected.xlsx", sheetName="EC", row.names=FALSE)
write.xlsx(sum_TD, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Final_RiskSNPOverlap_Corrected.xlsx", sheetName="TD", append=TRUE, row.names=FALSE)


#require(openxlsx)
#list_of_datasets <- list("EC_riskOverlaps" = sum_EC, "TD_riskOverlaps" = sum_TD)
#write.xlsx(list_of_datasets, file = "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/Final_RiskSNPOverlap.xlsx")



#################### UCSC format 


MakingUCSC <- function(TraitSum,vector,liftedCat,CNV_df){
  cat <- TraitSum[,c(1,5)]
  ecac <- as.data.frame(unique(vector))
  colnames(ecac) <- c("ID")
  for(i in 1:nrow(ecac)){
    ecac$Chr[i] <- unlist(strsplit(ecac$ID[i],split = "_"))[[1]]
    ecac$Start[i] <- unlist(strsplit(ecac$ID[i],split = "_"))[[2]]
    ecac$Stop[i] <- unlist(strsplit(ecac$ID[i],split = "_"))[[3]]
    ecac$Type[i] <- unlist(strsplit(ecac$ID[i],split = "_"))[[4]]
    ecac$Num.cases[i] <- length(unique(CNV_df$sampleid[CNV_df$cnv_id==ecac$ID[i] & CNV_df$cohort==1]))
    ecac$Num.controls[i] <- length(unique(CNV_df$sampleid[CNV_df$cnv_id==ecac$ID[i] & CNV_df$cohort==0]))
    ecac$Label[i] <- paste("Cases:",ecac$Num.cases[i],"_","Controls:",ecac$Num.controls[i],sep = "")
  }
  
  ecac_ucsc <- ecac[,c(2,3,4,8)]
  ecac_ucsc$Zero <- c("0")
  ecac_ucsc$Plus <- c("+")
  ecac_ucsc <- cbind(ecac_ucsc,ecac[,c(3:5)])
  ecac_ucsc$Chr <- paste("chr",ecac_ucsc$Chr,sep="")
  ecac_ucsc$Type <- ifelse(ecac_ucsc$Type=="del","255,0,0","0,0,255")
  
  
  SNP_ucsc <- as.data.frame(unique(TraitSum$RiskSNP))
  colnames(SNP_ucsc) <- c("Label")
  for(i in 1:nrow(SNP_ucsc)){
    SNP_ucsc$Chr[i] <- liftedCat$hg19.chr[liftedCat$hg19.rsID == SNP_ucsc$Label[i]]
    SNP_ucsc$Start[i] <- liftedCat$hg19.pos1[liftedCat$hg19.rsID == SNP_ucsc$Label[i]]
    SNP_ucsc$Stop[i] <- liftedCat$hg19.pos1[liftedCat$hg19.rsID == SNP_ucsc$Label[i]]
  }
  
  
  
  SNP_ucsc2 <- SNP_ucsc[,c(2,3,4,1)]
  SNP_ucsc2$Zero <- c("0")
  SNP_ucsc2$Plus <- c("+")
  SNP_ucsc2 <- cbind(SNP_ucsc2,SNP_ucsc[,c(3,4)])
  SNP_ucsc2$Type <- c("128,0,128")
  #SNP_ucsc2$Chr <- paste("chr",SNP_ucsc2$Chr,sep="")
  
  output <- rbind(SNP_ucsc2,ecac_ucsc)
  return(output)
}



t <- MakingUCSC(sum_EC,EC_vector,EC,ThreeProbeCNVs)


EC_UCSCFormat <- MakingUCSC(sum_EC,EC_vector,EC,ThreeProbeCNVs)
TD_UCSCFormat <- MakingUCSC(sum_TD,TD_vector,TD,ThreeProbeCNVs)

write.table(EC_UCSCFormat,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/ECoverlaps_redo.txt",quote = FALSE,row.names = FALSE,sep = '\t',col.names = FALSE)
write.table(TD_UCSCFormat,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_catalog_riskSNPs/TDoverlaps_redo.txt",quote = FALSE,row.names = FALSE,sep = '\t',col.names = FALSE)



