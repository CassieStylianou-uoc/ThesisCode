FreqMatch <- function(AssocOb,BinaryMatrix,SampleOb,case_n,control_n){
  for(i in 1:nrow(AssocOb)){
    Gene_data <- as.data.frame(BinaryMatrix[,which(colnames(BinaryMatrix) == AssocOb$name[i])])
    rownames(Gene_data) <- rownames(BinaryMatrix)
    Gene_cohort <- cbind(Gene_data,SampleOb$case_control[match(rownames(Gene_data),SampleOb$OncID)])
    colnames(Gene_cohort) <- c("CNV.presence","Cohort")
    cohortTable <- table(Gene_cohort$CNV.presence,Gene_cohort$Cohort)
    AssocOb$All.overlaps[i] <- sum(Gene_data)
    AssocOb$Case.overlaps[i] <- cohortTable[2,2]
    AssocOb$Case.prop[i] <- AssocOb$Case.overlaps[i]/case_n *100
    AssocOb$Control.overlaps[i] <- cohortTable[2,1]
    AssocOb$Control.prop[i] <- AssocOb$Control.overlaps[i]/control_n *100
    AssocOb$Gene.name[i] <- unlist(strsplit(AssocOb$name[i],"_"))[4]
    
    matrix_data <- as.data.frame(BinaryMatrix[,which(colnames(BinaryMatrix) == AssocOb$name[i])])
    rownames(matrix_data) <- rownames(BinaryMatrix)
    matrix_data_cohort <- cbind(matrix_data,SampleOb$case_control[match(rownames(Gene_data),SampleOb$OncID)],SampleOb$OncID[match(rownames(Gene_data),SampleOb$OncID)])
    colnames(matrix_data_cohort) <- c("CNV.presence","Cohort","sampleid")
    CNV_samples <- matrix_data_cohort[matrix_data_cohort$CNV.presence == 1,] #isolating the samples which have a CNV within gene
    case_samples <- paste(toString(rownames(CNV_samples[CNV_samples$Cohort == 1,])),sep = ";")
    control_samples <- paste(toString(rownames(CNV_samples[CNV_samples$Cohort == 0,])),sep = ";")
    AssocOb$Case.samples[i] <- case_samples
    AssocOb$Control.samples[i] <- control_samples
  }
  return(AssocOb)
  
}



binDir <- c("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/BinaryMatrixes/")
################ Reading in binary matrixes
b_2_del <- readRDS(paste0(binDir,"TwoProbe_del.rds"))
b_2_dup <- readRDS(paste0(binDir,"TwoProbe_dup.rds"))
b_2_lof <- readRDS(paste0(binDir,"TwoProbe_lof.rds"))

b_3_del <- readRDS(paste0(binDir,"ThreeProbe_del.rds"))
b_3_dup <- readRDS(paste0(binDir,"ThreeProbe_dup.rds"))
b_3_lof <- readRDS(paste0(binDir,"ThreeProbe_lof.rds"))

b_5_del <- readRDS(paste0(binDir,"FiveProbe_del.rds"))
b_5_dup <- readRDS(paste0(binDir,"FiveProbe_dup.rds"))
b_5_lof <- readRDS(paste0(binDir,"FiveProbe_lof.rds"))

b_10_del <- readRDS(paste0(binDir,"TenProbe_del.rds"))
b_10_dup <- readRDS(paste0(binDir,"TenProbe_dup.rds"))
b_10_lof <- readRDS(paste0(binDir,"TenProbe_lof.rds"))


assDir <- c("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/Associations/")
########### Reading in associations 

a_2_del <- readRDS(paste0(assDir,"TwoProbe_del_association.rds"))
a_2_dup <- readRDS(paste0(assDir,"TwoProbe_dup_association.rds"))
a_2_lof <- readRDS(paste0(assDir,"TwoProbe_lof_association.rds"))

a_3_del <- readRDS(paste0(assDir,"ThreeProbe_del_association.rds"))
a_3_dup <- readRDS(paste0(assDir,"ThreeProbe_dup_association.rds"))
a_3_lof <- readRDS(paste0(assDir,"ThreeProbe_lof_association.rds"))

a_5_del <- readRDS(paste0(assDir,"FiveProbe_del_association.rds"))
a_5_dup <- readRDS(paste0(assDir,"FiveProbe_dup_association.rds"))
a_5_lof <- readRDS(paste0(assDir,"FiveProbe_lof_association.rds"))

a_10_del <- readRDS(paste0(assDir,"TenProbe_del_association.rds"))
a_10_dup <- readRDS(paste0(assDir,"TenProbe_dup_association.rds"))
a_10_lof <- readRDS(paste0(assDir,"TenProbe_lof_association.rds"))



###################### Formatting all associations

GWAS_2_del <- FreqMatch(a_2_del,b_2_del,s,4115,17818)
GWAS_2_dup <- FreqMatch(a_2_dup,b_2_dup,s,4115,17818)
GWAS_2_lof <- FreqMatch(a_2_lof,b_2_lof,s,4115,17818)


GWAS_3_del <- FreqMatch(a_3_del,b_3_del,s,4115,17818)
GWAS_3_dup <- FreqMatch(a_3_dup,b_3_dup,s,4115,17818)
GWAS_3_lof <- FreqMatch(a_3_lof,b_3_lof,s,4115,17818)


GWAS_5_del <- FreqMatch(a_5_del,b_5_del,s,4115,17818)
GWAS_5_dup <- FreqMatch(a_5_dup,b_5_dup,s,4115,17818)
GWAS_5_lof <- FreqMatch(a_5_lof,b_5_lof,s,4115,17818)

GWAS_10_del <- FreqMatch(a_10_del,b_10_del,s,4115,17818)
GWAS_10_dup <- FreqMatch(a_10_dup,b_10_dup,s,4115,17818)
GWAS_10_lof <- FreqMatch(a_10_lof,b_10_lof,s,4115,17818)




GWAS_2_del <- GWAS_2_del[order(GWAS_2_del$p),]
GWAS_3_del <- GWAS_3_del[order(GWAS_3_del$p),]
GWAS_5_del <- GWAS_5_del[order(GWAS_5_del$p),]
GWAS_10_del <- GWAS_10_del[order(GWAS_10_del$p),]

GWAS_2_dup <- GWAS_2_dup[order(GWAS_2_dup$p),]
GWAS_3_dup <- GWAS_3_dup[order(GWAS_3_dup$p),]
GWAS_5_dup <- GWAS_5_dup[order(GWAS_5_dup$p),]
GWAS_10_dup <- GWAS_10_dup[order(GWAS_10_dup$p),]

GWAS_2_lof <- GWAS_2_lof[order(GWAS_2_lof$p),]
GWAS_3_lof <- GWAS_3_lof[order(GWAS_3_lof$p),]
GWAS_5_lof <- GWAS_5_lof[order(GWAS_5_lof$p),]
GWAS_10_lof <- GWAS_10_lof[order(GWAS_10_lof$p),]


####### Correcting - bonferroni



GWAS_2_del$p.adjust <- p.adjust(GWAS_2_del$p,"bonferroni")
GWAS_3_del$p.adjust <- p.adjust(GWAS_3_del$p,"bonferroni")
GWAS_5_del$p.adjust <- p.adjust(GWAS_5_del$p,"bonferroni")
GWAS_10_del$p.adjust <- p.adjust(GWAS_10_del$p,"bonferroni")


GWAS_2_dup$p.adjust <- p.adjust(GWAS_2_dup$p,"bonferroni")
GWAS_3_dup$p.adjust <- p.adjust(GWAS_3_dup$p,"bonferroni")
GWAS_5_dup$p.adjust <- p.adjust(GWAS_5_dup$p,"bonferroni")
GWAS_10_dup$p.adjust <- p.adjust(GWAS_10_dup$p,"bonferroni")

GWAS_2_lof$p.adjust <- p.adjust(GWAS_2_lof$p,"bonferroni")
GWAS_3_lof$p.adjust <- p.adjust(GWAS_3_lof$p,"bonferroni")
GWAS_5_lof$p.adjust <- p.adjust(GWAS_5_lof$p,"bonferroni")
GWAS_10_lof$p.adjust <- p.adjust(GWAS_10_lof$p,"bonferroni")

