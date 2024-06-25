### This script tests the sensitivity of the burden results to prove coverage threshold (i.e., across CNV sets) and to global CNV load limits


############ TWO PROBES

SS_2 <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SampleSummary_2probes.csv")
SS_2 <- SS_2[,-c(1)]
## How many samples had no CNVs (covered by three probes) - split by cohort
table(SS_2$case_control[SS_2$CNV==0]) #640 controls and 29 cases 

burdenI <- colnames(SS_2[c(3:23)])


burden_2 <- as.data.frame(burdenI)
for(i in 1:nrow(burden_2)){
  t <- t.test(SS_2[SS_2$case_control == 1,burdenI[i]],SS_2[SS_2$case_control == 0,burdenI[i]],alternative = "two.sided",var.equal = F)
  burden_2$Case.mean[i] <- t$estimate[1]
  burden_2$Case.sd[i] <- sd(SS_2[SS_2$case_control == 1,burdenI[i]],na.rm = FALSE)
  burden_2$Control.mean[i] <- t$estimate[2]
  burden_2$Control.sd[i] <- sd(SS_2[SS_2$case_control == 0,burdenI[i]],na.rm = FALSE)
  burden_2$Mean.difference[i] <- t$estimate[1] - t$estimate[2]
  burden_2$Lower.CI[i] <- t$conf.int[1]
  burden_2$Upper.CI[i] <- t$conf.int[2]
  burden_2$CI[i] <-  paste(t$conf.int[1],"to",t$conf.int[2])
  burden_2$P.value.full[i] <- t$p.value
  burden_2$P.value.refined[i] <- signif(t$p.value,digits = 3)
  print(i)
}

for(i in 1:nrow(burden_2)){
  burden_2$case_minus_control[i] <- burden_2$Case.mean[i] - burden_2$Control.mean[i]
  burden_2$case_control_ratio[i] <- burden_2$Case.mean[i] / burden_2$Control.mean[i]
}

############ FIVE PROBES

SS_5 <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SampleSummary_5probes.csv")
SS_5 <- SS_5[,-c(1)]
## How many samples had no CNVs (covered by three probes) - split by cohort
table(SS_5$case_control[SS_5$CNV==0]) #1468 controls and 195 cases 

burdenI <- colnames(SS_5[c(3:23)])


burden_5 <- as.data.frame(burdenI)
for(i in 1:nrow(burden_5)){
  t <- t.test(SS_5[SS_5$case_control == 1,burdenI[i]],SS_5[SS_5$case_control == 0,burdenI[i]],alternative = "two.sided",var.equal = F)
  burden_5$Case.mean[i] <- t$estimate[1]
  burden_5$Case.sd[i] <- sd(SS_5[SS_5$case_control == 1,burdenI[i]],na.rm = FALSE)
  burden_5$Control.mean[i] <- t$estimate[2]
  burden_5$Control.sd[i] <- sd(SS_5[SS_5$case_control == 0,burdenI[i]],na.rm = FALSE)
  burden_5$Mean.difference[i] <- t$estimate[1] - t$estimate[2]
  burden_5$Lower.CI[i] <- t$conf.int[1]
  burden_5$Upper.CI[i] <- t$conf.int[2]
  burden_5$CI[i] <-  paste(t$conf.int[1],"to",t$conf.int[2])
  burden_5$P.value.full[i] <- t$p.value
  burden_5$P.value.refined[i] <- signif(t$p.value,digits = 3)
  print(i)
}

for(i in 1:nrow(burden_5)){
  burden_5$case_minus_control[i] <- burden_5$Case.mean[i] - burden_5$Control.mean[i]
  burden_5$case_control_ratio[i] <- burden_5$Case.mean[i] / burden_5$Control.mean[i]
}


############ TEN PROBES

SS_10 <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SampleSummary_10probes.csv")
SS_10 <- SS_10[,-c(1)]
## How many samples had no CNVs (covered by three probes) - split by cohort
table(SS_10$case_control[SS_10$CNV==0]) #4090 controls and 771 cases 

burdenI <- colnames(SS_10[c(3:23)])


burden_10 <- as.data.frame(burdenI)
for(i in 1:nrow(burden_10)){
  t <- t.test(SS_10[SS_10$case_control == 1,burdenI[i]],SS_10[SS_10$case_control == 0,burdenI[i]],alternative = "two.sided",var.equal = F)
  burden_10$Case.mean[i] <- t$estimate[1]
  burden_10$Case.sd[i] <- sd(SS_10[SS_10$case_control == 1,burdenI[i]],na.rm = FALSE)
  burden_10$Control.mean[i] <- t$estimate[2]
  burden_10$Control.sd[i] <- sd(SS_10[SS_10$case_control == 0,burdenI[i]],na.rm = FALSE)
  burden_10$Mean.difference[i] <- t$estimate[1] - t$estimate[2]
  burden_10$Lower.CI[i] <- t$conf.int[1]
  burden_10$Upper.CI[i] <- t$conf.int[2]
  burden_10$CI[i] <-  paste(t$conf.int[1],"to",t$conf.int[2])
  burden_10$P.value.full[i] <- t$p.value
  burden_10$P.value.refined[i] <- signif(t$p.value,digits = 3)
  print(i)
}

for(i in 1:nrow(burden_10)){
  burden_10$case_minus_control[i] <- burden_10$Case.mean[i] - burden_10$Control.mean[i]
  burden_10$case_control_ratio[i] <- burden_10$Case.mean[i] / burden_10$Control.mean[i]
}



############ Write them all out to one excel file 

library(xlsx)
write.xlsx(burden_2, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SensitivityAnalysis.xlsx", sheetName="2_probes", row.names=FALSE)
write.xlsx(burden_5, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SensitivityAnalysis.xlsx", sheetName="5_probes", append=TRUE, row.names=FALSE)
write.xlsx(burden_10, file="C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SensitivityAnalysis.xlsx", sheetName="10_probes", append=TRUE, row.names=FALSE)





########################
### Reading in 3+ probe sample summary set 
SS <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/SampleSummary_3probe.rds")
table(SS$case_control)
burdenI <- colnames(SS[c(3:5)])

burdenFun <- function(SS_thresh,OriginalSS,burdenI){
  #### Excluding samples surpassing threshold
  sampleSummary <- OriginalSS[OriginalSS$CNV <= SS_thresh,]
  burdenOb <- as.data.frame(burdenI)
  
  ### performing burden
  for(i in 1:nrow(burdenOb)){
    t <- t.test(sampleSummary[sampleSummary$case_control == 1,burdenI[i]],sampleSummary[sampleSummary$case_control == 0,burdenI[i]],alternative = "two.sided",var.equal = F)
    burdenOb$Case.mean[i] <- t$estimate[1]
    burdenOb$Case.sd[i] <- sd(sampleSummary[sampleSummary$case_control == 1,burdenI[i]],na.rm = FALSE)
    burdenOb$Control.mean[i] <- t$estimate[2]
    burdenOb$Control.sd[i] <- sd(sampleSummary[sampleSummary$case_control == 0,burdenI[i]],na.rm = FALSE)
    burdenOb$Lower.CI[i] <- t$conf.int[1]
    burdenOb$Upper.CI[i] <- t$conf.int[2]
    burdenOb$CI[i] <- paste(t$conf.int[1]," to ",t$conf.int[2])
    burdenOb$P.value.full[i] <- t$p.value
    burdenOb$P.value.refined[i] <- signif(t$p.value,digits = 3)
    burdenOb$case_minus_control[i] <- burdenOb$Case.mean[i] - burdenOb$Control.mean[i]
    burdenOb$case_control_ratio[i] <- burdenOb$Case.mean[i] / burdenOb$Control.mean[i]
  }
  return(burdenOb)
}

Burden_50 <- burdenFun(50,SS,burdenI)
Burden_45 <- burdenFun(45,SS,burdenI)
Burden_40 <- burdenFun(40,SS,burdenI)
Burden_35 <- burdenFun(35,SS,burdenI)
Burden_30 <- burdenFun(30,SS,burdenI)
Burden_25 <- burdenFun(25,SS,burdenI)
Burden_20 <- burdenFun(20,SS,burdenI)
Burden_15 <- burdenFun(15,SS,burdenI)
Burden_10 <- burdenFun(10,SS,burdenI)




AdditionalFunc <- function(SS_thresh,OriginalSS){
  x <- OriginalSS[OriginalSS$CNV <= SS_thresh,]
  CaseExclude <- length(x$OncID[x$case_control == 1])
  ControlExclude <- length(x$OncID[x$case_control == 0])
  TotalSamplesExcluded <- length(x$OncID)
  print(paste("Cases excluded:",CaseExclude))
  print(paste("Controls excluded:",ControlExclude))
  print(paste("Total samples:",TotalSamplesExcluded))
}

AdditionalFunc(50,SS)
AdditionalFunc(45,SS)
AdditionalFunc(40,SS)
AdditionalFunc(35,SS)
AdditionalFunc(30,SS)
AdditionalFunc(25,SS)
AdditionalFunc(20,SS)
AdditionalFunc(15,SS)
AdditionalFunc(10,SS)



