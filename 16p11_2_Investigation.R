


##################################
## TESTS TO SEE IF BMI IS SIGNIFICANTLY DIFFERENT BETWEEN CARRIERS OF 16p11.2 DELETION AND NOT 


COI <- c("16_29595483_30156963_del")
c <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")
c_3 <- c[c$probes>=3,]
c_3$CNV_ID <- paste(c_3$chr,c_3$startpos,c_3$endpos,c_3$cnv_call,sep = "_")  

p16_c <- c_3[c_3$CNV_ID == COI,]

S <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_SampleSet.rds")

S_16p <- S[S$OncID %in% p16_c$sampleid,]

mean(S_16p$bmi[S_16p$bmi>=1]) #33.5
mean(S$bmi[S$bmi >=1 & S$OncID %!in% S_16p$OncID]) #26.55

S_withBMI <- S[S$bmi >1,] #16181 samples (13218 controls and 2963 cases)
bmi_with_deletion <- S_16p$bmi[S_16p$bmi >=1]
bmi_without_deletion <- S$bmi[S$bmi >=1 & S$OncID %!in% S_16p$OncID]
# Perform t-test
t_result <- t.test(bmi_with_deletion, bmi_without_deletion)

#######################################################



'%!in%' <- function(x,y)!('%in%'(x,y))
exp_ob <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/ExpressionNormalEndo_141candidateGenes.csv")
dosage_ob <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/ExpressionCancer_141candidateGenes_dosageSensitivity.csv")


ProxDel <- c("SPN",
             "QPRT",
             "C16orf54",
             "ZG16",
             "KIF22",
             "MAZ",
             "PRRT2",
             "PAGR1",
             "MVP",
             "CDIPT",
             "SEZ6L2",
             "ASPHD1",
             "KCTD13",
             "TMEM219",
             "TAOK2",
             "HIRIP3",
             "INO80E",
             "DOC2A",
             "C16orf92",
             "ALDOA",
             "PPP4C",
             "TBX6",
             "YPEL3",
             "GDPD3",
             "MAPK3")


exp_16 <- exp_ob[exp_ob$Gene.name %in% ProxDel,]
dosage_16 <- dosage_ob[dosage_ob$gene %in% ProxDel,]


DS_16p <- dosage_ob$gene[dosage_ob$pvalue < 0.0001 & dosage_ob$gene %in% ProxDel]


plot_16p <- ggplot(data=exp_ob[exp_ob$Gene.name %in% ProxDel,], aes(x=reorder(Gene.name,-TPM), y=TPM,fill = Gene.name %in% DS_16p)) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("grey80", "steelblue")) + theme_classic()+
  scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90))


zoomtest <- ggplot(data=exp_ob[exp_ob$Gene.name %in% ProxDel,], aes(x=reorder(Gene.name,-TPM), y=TPM,fill = Gene.name %in% DS_16p)) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("grey80", "steelblue")) + theme_classic()+
  scale_y_continuous(expand = c(0,0)) + theme(axis.text.x = element_text(angle = 90)) + coord_cartesian(xlim = c(19.2,25),ylim = c(0,10))



png("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/16p/Prox16p_barPlot.png",width = 700,height = 450)
plot_16p
dev.off()


png("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/16p/Prox16p_zoomedBar.png",width = 400,height = 400)
zoomtest
dev.off()


