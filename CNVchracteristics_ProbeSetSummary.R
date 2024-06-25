

### Reading in probe sets - labelled according to minimum probe length threshold 

two <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
three <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_3_200.rds")
five <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_5_200.rds")
ten <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_10_200.rds")

### Define a function to add CNV ID to all CNV data - this will allow for identification of unique CNVs 
formatFun <- function(cnv_PS){
  cnv_PS$chr <- paste("chr",cnv_PS$chr,sep="")
  cnv_PS$cnv_id <- paste(cnv_PS$chr,cnv_PS$startpos,cnv_PS$endpos,cnv_PS$cnv_call,sep="_")
  return(cnv_PS)
}


two <- formatFun(two)
three <- formatFun(three)
five <- formatFun(five)
ten <- formatFun(ten)


# PS_summary function definition
PS_summary <- function(cnv_ps) {
  # a: Number of deletions and duplications
  a <- cbind(nrow(cnv_ps[cnv_ps$cnv_call == "del", ]), nrow(cnv_ps[cnv_ps$cnv_call == "dup", ]))
  
  # b: Number of unique CNV IDs for deletions and duplications
  b <- cbind(length(unique(cnv_ps$cnv_id[cnv_ps$cnv_call == "del"])), length(unique(cnv_ps$cnv_id[cnv_ps$cnv_call == "dup"])))
  
  # c: Mean segment mean for deletions and duplications
  c <- cbind(round(mean(cnv_ps$segmean[cnv_ps$cnv_call == "del"]), 2), round(mean(cnv_ps$segmean[cnv_ps$cnv_call == "dup"]), 2))
  
  # dd: Median segment mean for deletions and duplications
  dd <- cbind(round(median(cnv_ps$segmean[cnv_ps$cnv_call == "del"]), 2), round(median(cnv_ps$segmean[cnv_ps$cnv_call == "dup"]), 2))
  
  # z: Range of segment means for deletions and duplications
  z <- cbind(paste(unlist(range(cnv_ps$segmean[cnv_ps$cnv_call == "del"]))[1], "to", unlist(range(cnv_ps$segmean[cnv_ps$cnv_call == "del"]))[2]), 
             paste(unlist(range(cnv_ps$segmean[cnv_ps$cnv_call == "dup"]))[1], "to", unlist(range(cnv_ps$segmean[cnv_ps$cnv_call == "dup"]))[2]))
  
  # d: Mean length of deletions and duplications (in kilobases)
  d <- cbind(round(mean(cnv_ps$length[cnv_ps$cnv_call == "del"]) / 1000, 2), round(mean(cnv_ps$length[cnv_ps$cnv_call == "dup"]) / 1000, 2))
  
  # e: Median length of deletions and duplications (in kilobases)
  e <- cbind(median(cnv_ps$length[cnv_ps$cnv_call == "del"]) / 1000, median(cnv_ps$length[cnv_ps$cnv_call == "dup"]) / 1000)
  
  # tempDat_case: Frequency table of sample IDs for cases (cohort == 1)
  tempDat_case <- as.data.frame(table(cnv_ps$sampleid[cnv_ps$cohort == 1]))
  
  # tempDat_control: Frequency table of sample IDs for controls (cohort == 0)
  tempDat_control <- as.data.frame(table(cnv_ps$sampleid[cnv_ps$cohort == 0]))
  
  # g: Mean number of probes for deletions and duplications
  g <- cbind(round(mean(cnv_ps$probes[cnv_ps$cnv_call == "del"]), 2), round(mean(cnv_ps$probes[cnv_ps$cnv_call == "dup"]), 2))
  
  # h: Median number of probes for deletions and duplications
  h <- cbind(median(cnv_ps$probes[cnv_ps$cnv_call == "del"]), median(cnv_ps$probes[cnv_ps$cnv_call == "dup"]))
  
  # singleton_del: Frequency table of CNV IDs for deletions
  singleton_del <- as.data.frame(table(cnv_ps$cnv_id[cnv_ps$cnv_call == "del"]))
  
  # singleton_dup: Frequency table of CNV IDs for duplications
  singleton_dup <- as.data.frame(table(cnv_ps$cnv_id[cnv_ps$cnv_call == "dup"]))
  
  # i: Percentage of singleton deletions and duplications (occurring only once)
  i <- cbind(round(nrow(singleton_del[singleton_del$Freq == 1, ]) / nrow(singleton_del) * 100, 2), round(nrow(singleton_dup[singleton_dup$Freq == 1, ]) / nrow(singleton_dup) * 100, 2))
  
  # f: Range of frequencies for deletions and duplications
  f <- cbind(paste(unlist(range(singleton_del$Freq))[1], "to", unlist(range(singleton_del$Freq))[2]), paste(unlist(range(singleton_dup$Freq))[1], "to", unlist(range(singleton_dup$Freq))[2]))
  
  # Return all the calculated summary statistics combined as rows
  return(rbind(a, b, c, dd, z, d, e, g, h, i, f))
}

# Apply the PS_summary function to multiple data sets and combine the results into a data frame
PS_Table <- as.data.frame(cbind(PS_summary(two), PS_summary(three), PS_summary(five), PS_summary(ten)))