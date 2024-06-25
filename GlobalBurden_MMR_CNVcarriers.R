#### This script determines global differences in CNV burden (split by type) for samples identified as carrying a CNV overlapping a MMR gene.
#### Also restricts samples to only include those that are predicted to carry a CNV assigned either PVS1 or PVS1_Strong in Chapter 7. 




library(readxl)
library(GenomicRanges)
aa <- read_excel("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/Thesis_writing/ThesisByPublication/VariantClassification/Phase 1/MANE/CG_VCmanuscript_transcriptCoordinates.xlsx",sheet=4)
bb <- aa[aa$Gene %in% c("MSH2","MLH1","PMS2","MSH6"),]


x <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
x$chr <- paste("chr",x$chr,sep="")
x$cnv_id <- paste(x$chr,x$startpos,x$endpos,x$cnv_call,sep="_")


bb_gr <- makeGRangesFromDataFrame(bb,seqnames.field = "Chr",start.field = "Varsome start",end.field = 'Calculated stop',keep.extra.columns = T)
x_gr <- makeGRangesFromDataFrame(x,seqnames.field = "chr",start.field = "startpos",end.field = "endpos",keep.extra.columns = T)

mmrMerge <- mergeByOverlaps(bb_gr,x_gr)
mmrMerge_df <- as.data.frame(mmrMerge)


### Pull out MMR carry sample ids
MMR_carriers <- unique(mmrMerge_df$sampleid)

############### Summarising for all 21,933 samples
S <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_SampleSet.rds")
fullSum <- as.data.frame(unique(S$OncID))
colnames(fullSum) <- c("sampleid")
c <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_3_200.rds")
c$chr <-  paste("chr",c$chr,sep="")
c$cnv_id <- paste(c$chr,c$startpos,c$endpos,c$cnv_call,sep="_")
for(i in 1:nrow(fullSum)){
  fullSum$CNV.count[i] <- length(unique(c$cnv_id[c$sampleid == fullSum$sampleid[i]]))
  fullSum$Del.count[i] <- length(unique(c$cnv_id[c$sampleid == fullSum$sampleid[i] & c$cnv_call=="del"]))
  fullSum$Dup.count[i] <- length(unique(c$cnv_id[c$sampleid == fullSum$sampleid[i] & c$cnv_call=="dup"]))
  print(i)
}


############### T-test: tests for significanct difference in global CNV burden between carriers and non-carries 

'%!in%' <- function(x,y)!('%in%'(x,y))
MMR_CNV <- t.test(fullSum$CNV.count[fullSum$sampleid %in% MMR_carriers],fullSum$CNV.count[fullSum$sampleid %!in% MMR_carriers],alternative = "two.sided",var.equal = F)
MMR_Del <-t.test(fullSum$Del.count[fullSum$sampleid %in% MMR_carriers],fullSum$Del.count[fullSum$sampleid %!in% MMR_carriers],alternative = "two.sided",var.equal = F)
MMR_Dup <-t.test(fullSum$Dup.count[fullSum$sampleid %in% MMR_carriers],fullSum$Dup.count[fullSum$sampleid %!in% MMR_carriers],alternative = "two.sided",var.equal = F)



########## Case - control comparison between carriers 
### All CNVs
Inter_carrier_CNV <- t.test(fullSum$CNV.count[fullSum$sampleid %in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$CNV.count[fullSum$sampleid %in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Inter_carrier_CNV$estimate[1] #Case carrier: 9.7
Inter_carrier_CNV$estimate[2] #Control carrier: 9.1
Inter_carrier_CNV$estimate[1] / Inter_carrier_CNV$estimate[2] #Case/control ratio: 1.06
Inter_carrier_CNV$estimate[1] - Inter_carrier_CNV$estimate[2] #Mean difference: 0.56
paste(round(Inter_carrier_CNV$conf.int[1],2)," to ",round(Inter_carrier_CNV$conf.int[2],2),sep="") #-2.73 to 3.86
Inter_carrier_CNV$p.value




## Duplications only
Inter_carrier_Dup <-t.test(fullSum$Dup.count[fullSum$sampleid %in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$Dup.count[fullSum$sampleid %in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Inter_carrier_Dup$estimate[1] #Case carrier: 4.58
Inter_carrier_Dup$estimate[2] #Control carrier: 4.94
Inter_carrier_Dup$estimate[1] / Inter_carrier_Dup$estimate[2] #Case/control ratio: 0.93
Inter_carrier_Dup$estimate[1] - Inter_carrier_Dup$estimate[2] #Mean difference: -0.36
paste(round(Inter_carrier_Dup$conf.int[1],2)," to ",round(Inter_carrier_Dup$conf.int[2],2),sep="") #-3.26 to 2.53
Inter_carrier_Dup$p.value





########## Case - control comparison between non-carriers 
Inter_Noncarrier_CNV <-t.test(fullSum$CNV.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$CNV.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Inter_Noncarrier_CNV$estimate[1] #Case carrier: 5.95
Inter_Noncarrier_CNV$estimate[2] #Control carrier: 4.89
Inter_Noncarrier_CNV$estimate[1] / Inter_Noncarrier_CNV$estimate[2] #Case/control ratio: 1.22
Inter_Noncarrier_CNV$estimate[1] - Inter_Noncarrier_CNV$estimate[2] #Mean difference: 1.06
paste(round(Inter_Noncarrier_CNV$conf.int[1],2)," to ",round(Inter_Noncarrier_CNV$conf.int[2],2),sep="") #0.93 to 1.18
Inter_Noncarrier_CNV$p.value




Inter_Noncarrier_Del <-t.test(fullSum$Del.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$Del.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Inter_Noncarrier_Del$estimate[1] #Case carrier: 3.23
Inter_Noncarrier_Del$estimate[2] #Control carrier: 2.8
Inter_Noncarrier_Del$estimate[1] / Inter_Noncarrier_Del$estimate[2] #Case/control ratio: 1.15
Inter_Noncarrier_Del$estimate[1] - Inter_Noncarrier_Del$estimate[2] #Mean difference: 0.43
paste(round(Inter_Noncarrier_Del$conf.int[1],2)," to ",round(Inter_Noncarrier_Del$conf.int[2],2),sep="") #0.35 to 0.52
Inter_Noncarrier_Del$p.value


Inter_Noncarrier_Dup <-t.test(fullSum$Dup.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$Dup.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)

Inter_Noncarrier_Dup$estimate[1] #Case carrier: 2.71
Inter_Noncarrier_Dup$estimate[2] #Control carrier: 2.08
Inter_Noncarrier_Dup$estimate[1] / Inter_Noncarrier_Dup$estimate[2] #Case/control ratio: 1.3
Inter_Noncarrier_Dup$estimate[1] - Inter_Noncarrier_Dup$estimate[2] #Mean difference: 0.63
paste(round(Inter_Noncarrier_Dup$conf.int[1],2)," to ",round(Inter_Noncarrier_Dup$conf.int[2],2),sep="") #0.54 to 0.71
Inter_Noncarrier_Dup$p.value



##### Outputting burden results to test for cohort difference between carriers and non-carrier groups

for(i in 1:nrow(burden_3)){
  t <- t.test(SS_3[SS_3$case_control == 1,burdenI[i]],SS_3[SS_3$case_control == 0,burdenI[i]],alternative = "two.sided",var.equal = F)
  burden_3$Case.mean[i] <- t$estimate[1]
  burden_3$Case.sd[i] <- sd(SS_3[SS_3$case_control == 1,burdenI[i]],na.rm = FALSE)
  burden_3$Control.mean[i] <- t$estimate[2]
  burden_3$Control.sd[i] <- sd(SS_3[SS_3$case_control == 0,burdenI[i]],na.rm = FALSE)
  burden_3$Mean.difference[i] <- t$estimate[1] - t$estimate[2]
  burden_3$Lower.CI[i] <- t$conf.int[1]
  burden_3$Upper.CI[i] <- t$conf.int[2]
  burden_3$CI[i] <-  paste(t$conf.int[1],"to",t$conf.int[2])
  burden_3$P.value.full[i] <- t$p.value
  burden_3$P.value.refined[i] <- signif(t$p.value,digits = 3)
  print(i)
}

for(i in 1:nrow(burden_3)){
  burden_3$case_minus_control[i] <- burden_3$Case.mean[i] - burden_3$Control.mean[i]
  burden_3$case_control_ratio[i] <- burden_3$Case.mean[i] / burden_3$Control.mean[i]
}





################## Testing to see if the effect holds up if all carrier samples are removed
case_samples <- S$OncID[S$case_control==1]
control_samples <- S$OncID[S$case_control==0]

CNV_noCarrier <- t.test(fullSum$CNV.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$CNV.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Del_noCarrier <- t.test(fullSum$Del.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$Del.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)
Dup_noCarrier <- t.test(fullSum$Dup.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% case_samples],fullSum$Dup.count[fullSum$sampleid %!in% MMR_carriers & fullSum$sampleid %in% control_samples],alternative = "two.sided",var.equal = F)





########################################## CHECKING TO SEE IF PVS1 OR PVS1_Strong weighting makes a difference

MMR_annotation <- read_excel("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_5/VC_AnnotationDocument.xlsx",sheet =1)
wantedMMR <- MMR_annotation$cnv_id[MMR_annotation$GENE %in% c("MSH2","MLH1","PMS2","MSH6") & MMR_annotation$PVS1_weighting %in% c("PVS1",'PVS1_Strong')]


x <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
x$chr <- paste("chr",x$chr,sep="")
x$cnv_id <- paste(x$chr,x$startpos,x$endpos,x$cnv_call,sep="_")


### need to use 2+ probe dataset - to allow for two probe CNVs that were annotated 
c_path <- x[x$cnv_id %in% wantedMMR,]

'%!in%' <- function(x,y)!('%in%'(x,y))
PathMMRcarrierVSnonCarrier_CNV <- t.test(fullSum$CNV.count[fullSum$sampleid %in% c_path$sampleid],fullSum$CNV.count[fullSum$CNV.count %!in% c_path$sampleid],alternative = "two.sided",var.equal = F)
PathMMRcarrierVSnonCarrier_DEL <- t.test(fullSum$Del.count[fullSum$sampleid %in% c_path$sampleid],fullSum$Del.count[fullSum$CNV.count %!in% c_path$sampleid],alternative = "two.sided",var.equal = F)
PathMMRcarrierVSnonCarrier_DUP <- t.test(fullSum$Dup.count[fullSum$sampleid %in% c_path$sampleid],fullSum$Dup.count[fullSum$CNV.count %!in% c_path$sampleid],alternative = "two.sided",var.equal = F)



cnv_summary <- fullSum



## compare MMR-CNV carriers to all others 
'%!in%' <- function(x,y)!('%in%'(x,y))
mean(cnv_summary$count[cnv_summary$sampleid %in% MMR_carriers]) #9.49
mean(cnv_summary$count[cnv_summary$sampleid %!in% MMR_carriers]) #5.08

t.test(cnv_summary$count[cnv_summary$sampleid %in% MMR_carriers],cnv_summary$count[cnv_summary$sampleid %!in% MMR_carriers],alternative = "two.sided",var.equal = F)



