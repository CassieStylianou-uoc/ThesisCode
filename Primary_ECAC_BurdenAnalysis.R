
SS_3 <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SampleSummary_ForBurden_3probes.csv")
SS_3 <- SS_3[,-c(1)]
## How many samples had no CNVs (covered by three probes) - split by cohort
table(SS_3$case_control[SS_3$CNV==0]) #761 controls and 45 cases 


burdenI <- colnames(SS_3[c(3:23)])


burden_3 <- as.data.frame(burdenI)



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

#write.csv(burden_3,"C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/3probeBurden.csv")

