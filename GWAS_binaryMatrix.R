################ Creates two functions -
### MakeBinMat:For deletions and duplications separately
### MakeBinMar: For loss of function CNVs - first determines level of overlap for each gene for all duplications and then only inludes partial gene duplications in GWAS. 




library(GenomicRanges)
library(tibble)
library(dplyr)
library(tidyr)

CNV <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
CNV$CNV_ID <- paste(CNV$chr,CNV$startpos,CNV$endpos,sep = "_")  
CNV$chr <- paste("chr",CNV$chr,sep="")


GENES <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/Thesis_writing/ThesisByPublication/GWAS_paper/R_GWAS_redo/UpdatedGeneList.csv")
GENES <- GENES[,-c(1,3,7)]
colnames(GENES) <- c("Gene","Start","Stop","Chr","Gene_ID")

## Converting genes list to GRanges object 
GENES_gr <- makeGRangesFromDataFrame(GENES,seqnames.field = "Chr",start.field = "Start",end.field = "Stop",keep.extra.columns = T,ignore.strand = T)
Samples <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_SampleSet.rds")

'%!in%' <- function(x,y)!('%in%'(x,y))


#### CNVobject = CNV
## GeneList_gr = GENES_gr
## SampleOb = Samples
MakeBinMat <- function(CNVobject,probeThresh,cnvType,GeneList_gr,SampleOb){
  cnv <- CNVobject[CNVobject$probes>= probeThresh & CNVobject$cnv_call == cnvType,]
  cnv_gr <- makeGRangesFromDataFrame(cnv,seqnames.field = "chr",start.field = "startpos",end.field = "endpos",keep.extra.columns = T,ignore.strand = T)
  merge <- mergeByOverlaps(GeneList_gr,cnv_gr)
  merge_df <- data.frame(sampleid = merge$sampleid,geneId = merge$Gene_ID,cnv_id = merge$CNV_ID)
  
  GeneSymbols <- merge_df$geneId[!duplicated(merge_df$geneId)]
  CNVids <- merge_df$cnv_id[!duplicated(merge_df$cnv_id)]
  
  binary_matrix <- merge_df %>%
    select(sampleid, geneId) %>%
    distinct() %>%
    mutate(value = 1) %>%
    spread(key = geneId, value = value, fill = 0)
  
  toAdd <- data.frame(SampleOb$OncID[SampleOb$OncID %!in% binary_matrix$sampleid])
  colnames(toAdd) <- c("sampleid")
  ZeroAdd <- matrix(0,ncol = ncol(binary_matrix)-1,nrow = nrow(toAdd))
  
  readytoAdd <- cbind(toAdd,ZeroAdd)
  colnames(readytoAdd)[2:ncol(readytoAdd)] <- colnames(binary_matrix)[2:ncol(binary_matrix)]
  binary_matrix2 <- rbind(binary_matrix,readytoAdd)
  binary_matrix2 <- column_to_rownames(binary_matrix2,var="sampleid")
  return(binary_matrix2)
}


testingBinMat <- MakeBinMat(CNV,3,"del",GENES_gr,Samples)


########### Updating function to do loss of function things

GeneList_LOF <- GENES
for(i in 1:nrow(GeneList_LOF)){
  GeneList_LOF$Gene.id[i] <- paste(unlist(strsplit(GeneList_LOF$Gene_ID[i],"_"))[1:3],collapse="_") 
}

GeneList_LOF_gr <- makeGRangesFromDataFrame(GeneList_LOF,seqnames.field = "Chr",start.field = "Start",end.field = "Stop",keep.extra.columns = T,ignore.strand = T)

MakeBinMat_LOF <- function(CNVobject,probeThresh,GeneList_gr,SampleOb){
  cnvs <- CNVobject[CNVobject$probes>= probeThresh,]
  cnv_gr <- makeGRangesFromDataFrame(cnvs,seqnames.field = "chr",start.field = "startpos",end.field = "endpos",keep.extra.columns = T,ignore.strand = T)
  ## Assessing overlap between CNVs and gene list 
  t <- olRanges(query = cnv_gr,subject = GeneList_gr,output = "gr")
  t_df <- data.frame(sampleid = t$sampleid,geneChr = t$space,geneStart = t$Sstart,geneStop = t$Send,cnvid = t$CNV_ID,overlapType = t$OLtype,cnvType = t$cnv_call)
  t_df$Gene.id <- paste(t_df$geneChr,t_df$geneStart,t_df$geneStop,sep="_")
  t_full <- merge(t_df,GeneList_LOF,by = "Gene.id")
  t_full <- t_full[,-c(10,11,12)]
  wantedCNVs <- t_full[t_full$overlapType %in% c("olup","contained","oldown") & t_full$cnvType == "dup" | t_full$cnvType=="del",]
  
  
  GeneSymbols <- unique(wantedCNVs$Gene_ID)
  CNVids <- unique(wantedCNVs$cnvid)
  
  # binary_matrix <- t_full %>%
  #   select(sampleid,Gene_ID) %>%
  #   distinct() %>%
  #   mutate(value = 1) %>%
  #   spread(key = Gene_ID,value = value,fill = 0)
  
  
  binary_matrix <- wantedCNVs %>%
    select(sampleid,Gene_ID) %>%
    distinct() %>%
    mutate(value = 1) %>%
    spread(key = Gene_ID,value = value,fill = 0)
  
  
  toAdd <- data.frame(SampleOb$OncID[SampleOb$OncID %!in% binary_matrix$sampleid])
  colnames(toAdd) <- c("sampleid")
  ZeroAdd <- matrix(0,ncol = ncol(binary_matrix)-1,nrow = nrow(toAdd))
  readytoAdd <- cbind(toAdd,ZeroAdd)
  colnames(readytoAdd)[2:ncol(readytoAdd)] <- colnames(binary_matrix)[2:ncol(binary_matrix)]
  binary_matrix2 <- rbind(binary_matrix,readytoAdd)
  binary_matrix2 <- column_to_rownames(binary_matrix2,var="sampleid")
  return(binary_matrix2)
  
}


################### Creating binary matrix's for GWAS
BinMatDir <- c("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/BinaryMatrixes/")

## TWO PROBE
TwoProbe_del <-  MakeBinMat(CNV,2,"del",GENES_gr,Samples)
TwoProbe_dup <-  MakeBinMat(CNV,2,"dup",GENES_gr,Samples)
TwoProbe_lof <-  MakeBinMat_LOF(CNV,2,GeneList_LOF_gr,Samples)

saveRDS(TwoProbe_del,paste0(BinMatDir,"TwoProbe_del.rds"))
saveRDS(TwoProbe_dup,paste0(BinMatDir,"TwoProbe_dup.rds"))
saveRDS(TwoProbe_lof,paste0(BinMatDir,"TwoProbe_lof.rds"))


## THREE PROBE
ThreeProbe_del <-  MakeBinMat(CNV,3,"del",GENES_gr,Samples)
ThreeProbe_dup <-  MakeBinMat(CNV,3,"dup",GENES_gr,Samples)
ThreeProbe_lof <-  MakeBinMat_LOF(CNV,3,GeneList_LOF_gr,Samples)

saveRDS(ThreeProbe_del,paste0(BinMatDir,"ThreeProbe_del.rds"))
saveRDS(ThreeProbe_dup,paste0(BinMatDir,"ThreeProbe_dup.rds"))
saveRDS(ThreeProbe_lof,paste0(BinMatDir,"ThreeProbe_lof.rds"))


## FIVE PROBE
FiveProbe_del <-  MakeBinMat(CNV,5,"del",GENES_gr,Samples)
FiveProbe_dup <-  MakeBinMat(CNV,5,"dup",GENES_gr,Samples)
FiveProbe_lof <-  MakeBinMat_LOF(CNV,5,GeneList_LOF_gr,Samples)


saveRDS(FiveProbe_del,paste0(BinMatDir,"FiveProbe_del.rds"))
saveRDS(FiveProbe_dup,paste0(BinMatDir,"FiveProbe_dup.rds"))
saveRDS(FiveProbe_lof,paste0(BinMatDir,"FiveProbe_lof.rds"))



## TEN PROBE
TenProbe_del <-  MakeBinMat(CNV,10,"del",GENES_gr,Samples)
TenProbe_dup <-  MakeBinMat(CNV,10,"dup",GENES_gr,Samples)
TenProbe_lof <-  MakeBinMat_LOF(CNV,10,GeneList_LOF_gr,Samples)


saveRDS(TenProbe_del,paste0(BinMatDir,"TenProbe_del.rds"))
saveRDS(TenProbe_dup,paste0(BinMatDir,"TenProbe_dup.rds"))
saveRDS(TenProbe_lof,paste0(BinMatDir,"TenProbe_lof.rds"))









pls_xxx <- column_to_rownames(xxx,var = "sampleid")




### Creating function that means partial duplication overlap can be determined
olRanges <- function(query, subject, output="gr", ...) {
  require(GenomicRanges); require(IRanges)
  
  ## Input check
  if(!((class(query)=="GRanges" & class(subject)=="GRanges") | (class(query)=="IRanges" & class(subject)=="IRanges"))) {
    stop("Query and subject need to be of same class, either GRanges or IRanges!")
  }
  
  ## Find overlapping ranges
  if(class(query)=="GRanges") {
    seqlengths(query) <- rep(NA, length(seqlengths(query)))
    seqlengths(subject) <- rep(NA, length(seqlengths(subject)))
  }
  olindex <- as.matrix(findOverlaps(query, subject, ...))
  query <- query[olindex[,1]]
  subject <- subject[olindex[,2]]
  olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
  
  ## Pre-queries for overlaps
  startup <- olma[,"Sstart"] < olma[,"Qstart"]
  enddown <- olma[,"Send"] > olma[,"Qend"]
  startin <- olma[,"Sstart"] >= olma[,"Qstart"] & olma[,"Sstart"] <= olma[,"Qend"]
  endin <- olma[,"Send"] >= olma[,"Qstart"] & olma[,"Send"] <=  olma[,"Qend"]
  
  ## Overlap types
  olup <- startup & endin
  oldown <- startin & enddown
  inside <- startin & endin 
  contained <- startup & enddown
  
  ## Overlap types in one vector
  OLtype <- rep("", length(olma[,"Qstart"]))
  OLtype[olup] <- "olup"
  OLtype[oldown] <- "oldown"
  OLtype[inside] <- "inside" 
  OLtype[contained] <- "contained"
  
  ## Overlap positions
  OLstart <- rep(0, length(olma[,"Qstart"]))
  OLend <- rep(0, length(olma[,"Qstart"]))
  OLstart[olup] <- olma[,"Qstart"][olup]
  OLend[olup] <- olma[,"Send"][olup]
  OLstart[oldown] <- olma[,"Sstart"][oldown]
  OLend[oldown] <- olma[,"Qend"][oldown]
  OLstart[inside] <- olma[,"Sstart"][inside]
  OLend[inside] <- olma[,"Send"][inside]
  OLstart[contained] <- olma[,"Qstart"][contained]
  OLend[contained] <- olma[,"Qend"][contained]
  
  ## Absolute and relative length of overlaps
  OLlength <- (OLend - OLstart) + 1
  OLpercQ <- OLlength/width(query)*100
  OLpercS <- OLlength/width(subject)*100
  
  ## Output type
  oldf <- data.frame(Qindex=olindex[,1], Sindex=olindex[,2], olma, OLstart, OLend, OLlength, OLpercQ, OLpercS, OLtype)
  if(class(query) == "GRanges") {
    oldf <- cbind(space=as.character(seqnames(query)), oldf)
  }
  if(output=="df") {
    return(oldf)
  }
  if(output=="gr") {
    if(class(query)=="GRanges") {
      elementMetadata(query) <- cbind(as.data.frame(elementMetadata(query)), oldf)
    }
    if(class(query)=="IRanges") {
      query <- GRanges(seqnames = Rle(rep("dummy", length(query))), ranges = IRanges(start=oldf[,"Qstart"], end=oldf[,"Qend"]), strand = Rle(strand(rep("+", length(query)))), oldf)  
    }
    return(query)
  }
}


