library(GenomicRanges)


############# OVERLAPPING RANGES FUNCTION
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






# Reading in ECAC CNVs
cnvs <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
cnvs$chr <- paste("chr",cnvs$chr,sep="")
cnvs$CNV_ID <- paste(cnvs$chr,cnvs$startpos,cnvs$endpos,cnvs$cnv_call,sep="_")
# Making a unique set of CNVs for annotation 
u_cnvs <- cnvs[!duplicated(cnvs$CNV_ID),c(3,4,5,6,8,10,16)]


##Setting directory to read in genomic context files (in BED format)
dir1 <- c("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/GeneExon_data/")
## Reading in gene list  - 17601 in total 
g <- readRDS(paste0(dir1,"2023_geneList.rds"))
## making gene ID - concatenating chromosome, start and stop position. 
g$Gene.id <- paste(g$chromosome_name,g$start_position,g$end_position,sep="_")

## Reading in exon data - exon info for 15,242 genes 
e <- readRDS(paste0(dir1,"2023_exonList.rds"))
e$chromosome_name <- paste("chr",e$chromosome_name,sep="")
e$chromosome_name <- gsub("chrX","chr23",e$chromosome_name)
e <- e[,-c(4)]

## Converting gene and exon data to a GenomicRanges to allow overlap assessment 
genes_gr <- makeGRangesFromDataFrame(g,keep.extra.columns = T,seqnames.field = "chromosome_name",start.field = "start_position",end.field = "end_position")
exons_gr <- makeGRangesFromDataFrame(e,keep.extra.columns = T, seqnames.field = "chromosome_name",start.field = "exon_chrom_start",end.field ="exon_chrom_end")
cnv_gr <- makeGRangesFromDataFrame(u_cnvs,keep.extra.columns = T, seqnames.field = "chr",start.field = "startpos",end.field = "endpos")

######### Performing merge function with gene list and all CNVs to use in-conjunction with annotation object
GeneMerge <- mergeByOverlaps(cnv_gr,genes_gr)

listX <- list()
for(i in 1:nrow(GeneMerge)){
  listX[[i]] <- as.data.frame(GeneMerge[i,])
}

GeneMerge_df <- do.call(rbind,listX)


############
#### Annotates each unique CNV with the the following - uses the olRanges function defined aove to determine the level of overlap between genes 
# Total number of genes overlapped
# Which genes are overlapped
# Total count of genes completely overlapped
# Which genes are completely overlapped
## Else condition:
# Assigns zero for all columns in that row if no gene overlap is predicted for any given CNV. 

#65,261
for(i in 1:nrow(u_cnvs)){
  tempO <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == u_cnvs$CNV_ID[i],],genes_gr,output="gr")
  if(length(tempO$OLtype)>=1){
    temp_df <- as.data.frame(tempO,row.names = NULL)
    temp_df$Gene.ID <- paste(temp_df$space,temp_df$Sstart,temp_df$Send,sep="_")
    for(j in 1:nrow(temp_df)){
      temp_df$Gene.name[j] <- g$hgnc_symbol[g$Gene.id == temp_df$Gene.ID[j]]
    }      
    u_cnvs$ID_confirm[i] <- temp_df$CNV_ID[1]
    u_cnvs$Gene.count[i] <- nrow(temp_df)
    u_cnvs$Genes.overlapped[i] <- paste(temp_df$Gene.name,collapse = ", ")
    u_cnvs$Whole.gene.count[i] <- nrow(temp_df[temp_df$OLpercS == 100,]) 
    u_cnvs$Genes.completelyOverlapped[i] <- paste(temp_df$Gene.name[temp_df$OLpercS == 100],collapse = ", ")
  }else{
    u_cnvs$ID_confirm[i] <- cnv_gr$CNV_ID[cnv_gr$CNV_ID == u_cnvs$CNV_ID[i]]
    u_cnvs$Gene.count[i] <- c("0")
    u_cnvs$Genes.overlapped[i] <- c("0")
    u_cnvs$Whole.gene.count[i] <- c("0")
    u_cnvs$Genes.completelyOverlapped[i] <- c("0")
  }
  print(i)
}

############
#### Annotates each unique CNV with the the following - uses the olRanges function defined aove to determine the level of overlap between exons 
# Total number of genes overlapped
# Which genes are overlapped
# Total count of genes completely overlapped
# Which genes are completely overlapped
## Else condition:
# Assigns zero for all columns in that row if no gene overlap is predicted for any given CNV. 


exon_cnvs <- cnvs[!duplicated(cnvs$CNV_ID),c(3,4,5,6,8,10,16)]
exonTable <- as.data.frame(table(e$Gene.name))
colnames(exonTable) <-c("Gene","Exon.count")

#65,261
nrow(exon_cnvs)
for(i in 1:nrow(exon_cnvs)){
  temp <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == exon_cnvs$CNV_ID[i],],exons_gr,output="gr")
  if(length(temp$OLtype) >=1){
    temp_df <- as.data.frame(temp,row.names=NULL)
    temp_df$Exon.id <- paste(temp_df$space,temp_df$Sstart,temp_df$Send,sep="_")
    for(j in 1:nrow(temp_df)){
      temp_df$Gene.name[j] <- e$Gene.name[e$exon_id == temp_df$Exon.id[j]]
    }
    exon_cnvs$Exonic.Gene.Count[i] <- length(unique(temp_df$Gene.name))
    exon_cnvs$Gene.overlapped.exonic[i] <- paste(unique(temp_df$Gene.name), collapse = ", ")
    exon_cnvs$Total.exon.count[i] <- nrow(temp_df[!duplicated(temp_df$Exon.id),])
    exon_cnvs$Whole.exon.count[i] <- nrow(temp_df[temp_df$OLpercS == 100,])
    
    tab <- as.data.frame(table(temp_df$Gene.name))
    colnames(tab) <-c("Gene","Exon.count.overlapped")
    merge <- merge(tab,exonTable,by = "Gene")
    
    exon_cnvs$Whole.gene.overlap.exonic[i] <- nrow(merge[merge$Exon.count.overlapped == merge$Exon.count,])
    exon_cnvs$Exonic.wholeOverlap[i] <- paste(merge$Gene[merge$Exon.count.overlapped == merge$Exon.count],collapse = ", ")
    exon_cnvs$Exon.Profile[i] <- paste(merge$Gene,":",merge$Exon.count.overlapped,"of",merge$Exon.count,collapse = " /// ") ## For each gene in which exonic overlap is predicted. Provides the ratio of how many exons are involved in the CNV and teh total number of exons in any given genes. 
  }else{
    exon_cnvs$Exonic.Gene.Count[i] <- c("0")
    exon_cnvs$Gene.overlapped.exonic[i] <- c("0")
    exon_cnvs$Total.exon.count[i]<- c("0")
    exon_cnvs$Whole.exon.count[i]<- c("0")
    exon_cnvs$Whole.gene.overlap.exonic[i] <- c("0")
    exon_cnvs$Exonic.wholeOverlap[i] <- c("0")
    exon_cnvs$Exon.Profile[i] <- c("0")
  }
  print(i)
}


## Annotetes CNV for overlap with a gene with a Pli score > 0.9 

for( i in 1:nrow(currentAnno)){
  anyOverlap <- unlist(strsplit(currentAnno$Genes.overlapped[i],split=", "))
  completeOverlap <- unlist(strsplit(currentAnno$Genes.completelyOverlapped[i],split=", "))
  
  currentAnno$pli_9_anyOverlap[i] <- length(intersect(anyOverlap,pli_0.9))
  currentAnno$pli_9_CompleteOverlap[i] <- length(intersect(completeOverlap,pli_0.9))
  
  currentAnno$pli_9_genesOverlapped[i] <- paste(intersect(anyOverlap,pli_0.9),collapse = ", ")
  currentAnno$pli_9_genesCompletelyOverlapped[i] <- paste(intersect(completeOverlap,pli_0.9),collapse=", ")
  
  print(i)
}

######### Adding frequency information to CNVs 
#65,261 CNVs 
for(i in 1:nrow(currentAnno)){
  Cdata <- cnvs[cnvs$CNV_ID == currentAnno$CNV_ID[i],]
  currentAnno$Freq.count[i] <- length(unique(Cdata$sampleid))
  currentAnno$Case.freq.count[i] <- length(unique(Cdata$sampleid[Cdata$cohort==1]))  
  currentAnno$Control.freq.count[i] <- length(unique(Cdata$sampleid[Cdata$cohort==0]))  
  print(i)
}




############# reading in CPG/DNaseHS/and miRNA data (hg19)

CpG <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_1/AnnotationData/RDS_GenomicContext_Files/CpG_gr.RDS")
DNase <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_1/AnnotationData/RDS_GenomicContext_Files/DNaseHS_gr.RDS")
miRNA <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_1/AnnotationData/RDS_GenomicContext_Files/mi_sno_gr.RDS")



library(GenomicRanges)

CpG_merge <- mergeByOverlaps(cnv_gr,CpG) #30,671
DNase_merge <- mergeByOverlaps(cnv_gr,DNase) #3235545
miRNA_merge <- mergeByOverlaps(cnv_gr,miRNA) #2,512


CPG_slot<- CpG_merge@listData
if (!inherits(CPG_slot, "data.frame")) {
  CPG_merge_df <- data.frame(CPG_slot)
}

#DNase_slot <- DNase_merge@listData
#if (!inherits(DNase_slot, "data.frame")) {
#  DNase_merge_df <- data.frame(DNase_slot)
#}

miRNA_slot <- miRNA_merge@listData
if (!inherits(miRNA_slot, "data.frame")) {
  miRNA_merge_df <- data.frame(miRNA_slot)
}

CPG_merge_df$CpGID <- paste(CPG_merge_df$CpG.seqnames,CPG_merge_df$CpG.start,CPG_merge_df$CpG.end,sep="_")
length(unique(CPG_merge_df$CpGID)) #12,349
DNase_merge_df$DnaseID <- paste(DNase_merge_df$DNase.seqnames,DNase_merge_df$DNase.start,DNase_merge_df$DNase.end,sep="_")
length(unique(DNase_merge_df$DnaseID)) #1,370,092
miRNA_merge_df$ID <- paste(miRNA_merge_df$miRNA.seqnames,miRNA_merge_df$miRNA.start,miRNA_merge_df$miRNA.end,miRNA_merge_df$type,sep="_")
length(unique(miRNA_merge_df$ID)) #1037

for(i in 1:nrow(currentAnno)){
  c <- CPG_merge_df[CPG_merge_df$CNV_ID == currentAnno$CNV_ID[i],]
  d <- DNase_merge_df[DNase_merge_df$CNV_ID == currentAnno$CNV_ID[i],]
  m <- miRNA_merge_df[miRNA_merge_df$CNV_ID == currentAnno$CNV_ID[i],]
  
  currentAnno$CpG.overlaps[i] <- length(unique(c$CpGID))
  currentAnno$DNase.overlaps[i] <- length(unique(d$DnaseID))
  currentAnno$miRNA.overlaps[i] <- length(unique(m$ID[m$type == "miRNA"]))
  currentAnno$CDBox.overlaps[i] <- length(unique(m$ID[m$type == "CDBox"]))
  currentAnno$HAcaBox.overlaps[i] <- length(unique(m$ID[m$type == "HAcaBox"]))
  currentAnno$scaRna.overlaps[i] <- length(unique(m$ID[m$type == "scaRna"]))
  
  temp_c <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],CpG,output = "gr")
  temp_d <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],DNase,output = "gr")
  temp_m <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],miRNA[miRNA$type=="miRNA"],output = "gr")
  temp_CDBox <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],miRNA[miRNA$type=="CDBox"],output = "gr")
  temp_HacaBox <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],miRNA[miRNA$type=="HAcaBox"],output = "gr")
  temp_scaRna <- olRanges(query = cnv_gr[cnv_gr$CNV_ID == currentAnno$CNV_ID[i],],miRNA[miRNA$type=="scaRna"],output = "gr")
  
  x<-temp_c@elementMetadata
  if (!inherits(x, "data.frame")) {
    x_df <- data.frame(x)
  }
  
  y <- temp_d@elementMetadata
  if (!inherits(y, "data.frame")) {
   y_df <- data.frame(y)
  }
  z <- temp_m@elementMetadata
  if (!inherits(z, "data.frame")) {
    z_df <- data.frame(z)
  }
  
  
  currentAnno$CpG.completeOverlaps[i] <- nrow(x[x$OLpercS == 100,])
  #currentAnno$DNase.completeOverlaps[i] <- nrow(y[y$OLpercS == 100,])
  currentAnno$miRNA.completeOverlaps[i] <- nrow(z[z$OLpercS == 100,])
  
  
  xx <- temp_CDBox@elementMetadata
  if(!inherits(xx,"data.frame")){
    xx_df <- data.frame(xx)
  }
  
  yy <- temp_HacaBox@elementMetadata
  if(!inherits(yy,"data.frame")){
    yy_df <- data.frame(yy)
  }
  
  zz <- temp_scaRna@elementMetadata
  if(!inherits(zz,"data.frame")){
    zz_df <- data.frame(zz)
  }
  
  currentAnno$CDBox.completeOverlaps[i] <- nrow(xx[xx$OLpercS == 100,])
  currentAnno$HAcaBox.completeOverlaps[i] <- nrow(yy[yy$OLpercS == 100,])
  currentAnno$scaRna.completeOverlaps[i] <- nrow(zz[zz$OLpercS == 100,])
  
  print(i)
}


saveRDS(currentAnno,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_1/ECAC_CNV_Annotation.rds")





