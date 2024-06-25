
library(GenomicRanges)
library(readxl)

g <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/Genes_final.rds")
g_gr <- makeGRangesFromDataFrame(g,keep.extra.columns = T,seqnames.field = "chromosome_name",start.field = "start_position",end.field = "end_position")

DNA_repair <- read_excel("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/HumanDNArepairGenes.xlsx",sheet = 2)

DNArepair_data <- g[g$hgnc_symbol %in% DNA_repair$DNA_repairGenes,]
#write.csv(DNArepair_data,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/HumanDNARepair_geneData.csv")
DNArepair_data <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/HumanDNARepair_geneData.csv")

nrow(DNArepair_data) #206 genes 

DNArepair_data_gr <- makeGRangesFromDataFrame(DNArepair_data,seqnames.field = "chromosome_name",keep.extra.columns = T,start.field = "start_position",end.field = "end_position")

cnvs <-  readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_3_200.rds")
cnvs$chr <- paste("chr",cnvs$chr,sep="")
cnvs$cnv_id <- paste(cnvs$chr,cnvs$startpos,cnvs$endpos,cnvs$cnv_call,sep="_")
cnvs_gr <- makeGRangesFromDataFrame(cnvs,seqnames.field = "chr",start.field = "startpos",end.field = "endpos",keep.extra.columns = T)
samples <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_SampleSet.rds")



###### Finding overlaps with CNVs and DNA repair genes

m <- mergeByOverlaps(DNArepair_data_gr,cnvs_gr)
list <- list()
for(i in 1:nrow(m)){
  list[[i]] <- as.data.frame(m[i,])
  print(i)
}

m_df <- do.call(rbind,list)


#######################

COI <- m_df

sampleTable <- as.data.frame(table(COI$sampleid))

for(i in 1:nrow(sampleTable)){
  d <- COI[COI$sampleid==sampleTable$Var1[i],]
  sampleTable$Gene[i] <- paste(d$hgnc_symbol,collapse = ", ")
  sampleTable$CNVs[i] <- paste(d$cnv_id,collapse = ", ")
  sampleTable$Cohort[i] <- d$cohort[1]
}

MMRcarrier_IDS <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/MismatchRepairGene_carriers.csv")
MMRcarrier_IDS <- MMRcarrier_IDS[,-c(1)]


### Pull out MMR carry sample ids
DNArepair_carriers <- unique(m_df$sampleid)


DNArepair_carriers <- unique(append(DNArepair_carriers,MMRcarrier_IDS))
length(DNArepair_carriers)


samples <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_SampleSet.rds")
sample_tally <- as.data.frame(samples$OncID)
colnames(sample_tally) <- c("OncID")
for(i in 1:nrow(sample_tally)){
  sample_tally$CNV.count[i] <- length(unique(cnvs$cnv_id[cnvs$sampleid == sample_tally$OncID[i]]))
  sample_tally$Del.count[i] <- length(unique(cnvs$cnv_id[cnvs$sampleid == sample_tally$OncID[i] & cnvs$cnv_call=="del"]))
  sample_tally$Dup.count[i] <- length(unique(cnvs$cnv_id[cnvs$sampleid == sample_tally$OncID[i] & cnvs$cnv_call=="dup"]))
}


'%!in%' <- function(x,y)!('%in%'(x,y))
DnaRepair_cnv <- t.test(sample_tally$CNV.count[sample_tally$OncID %in% DNArepair_carriers],sample_tally$CNV.count[sample_tally$OncID %!in% DNArepair_carriers],alternative = "two.sided",var.equal = F)
DnaRepair_del <-t.test(sample_tally$Del.count[sample_tally$OncID %in% DNArepair_carriers],sample_tally$Del.count[sample_tally$OncID %!in% DNArepair_carriers],alternative = "two.sided",var.equal = F)
DnaRepair_dup <-t.test(sample_tally$Dup.count[sample_tally$OncID %in% DNArepair_carriers],sample_tally$Dup.count[sample_tally$OncID %!in% DNArepair_carriers],alternative = "two.sided",var.equal = F)

case_samples <- samples$OncID[samples$case_control==1]
control_samples <- samples$OncID[samples$case_control==0]

########## Case - control comparison between carriers 
carrier_all <- t.test(sample_tally$CNV.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$CNV.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)
carrier_del <-t.test(sample_tally$Del.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Del.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)
carrier_dup <-t.test(sample_tally$Dup.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Dup.count[sample_tally$OncID %in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)


########## Case - control comparison between non-carriers 
non_carrier_cnv <- t.test(sample_tally$CNV.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$CNV.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)

non_carrier_del <- t.test(sample_tally$Del.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Del.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)

non_carrier_dup <- t.test(sample_tally$Dup.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Dup.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)



################## Testing to see if the effect holds up if all carrier samples are removed

noCarrier_CNV <- t.test(sample_tally$CNV.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$CNV.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)
noCarrier_Del <- t.test(sample_tally$Del.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Del.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)
noCarrier_Dup <- t.test(sample_tally$Dup.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% case_samples],sample_tally$Dup.count[sample_tally$OncID %!in% DNArepair_carriers & sample_tally$OncID %in% control_samples],alternative = "two.sided",var.equal = F)






