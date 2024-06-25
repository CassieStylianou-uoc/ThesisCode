
### libraries needed 
library(meta)
library(readxl)
library(dmetar)

## Reading in sheet from excel sheet with all data necessary on it 
meta <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/Meta_analysisData_2023.xlsx",sheet=1)

## Random effects model used:

meta_RE_SMD <- metacont(Ne,Me,Se,Nc,Mc,Sc, data = meta,studlab =paste(Author),fixed = FALSE,random = TRUE,prediction = TRUE,sm="SMD",method.smd = "Cohen")
## summary for meta 
summary(meta_RE_SMD)
forest(meta_RE_SMD,digits.se = 1,digits.mean = 1,sortvar = TE,leftcols = c("Cancer type","Author"),leftlabs = c("Cancer type","Author"))

m.gen.Influence <- InfluenceAnalysis(meta_RE_SMD,random = TRUE)

png(file = "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/InfluenceAnalysis_all.png", width = 3000, height = 2400, res = 300)
plot(m.gen.Influence,"es")
dev.off()


###### Funnel plot - all studies
png(file = "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/FunnelPlot_noLabels_all.png", width = 2500, height = 2000, res = 300)
col.contour = c("lavender", "lightblue1", "lightcyan1")
funnel.meta(meta_RE_SMD, xlim = c(-0.6, 1),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour,studlab = FALSE)
legend(x = 0.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)
dev.off()
#####################




## Meta analysis - excluding both dataset from Park
meta_noPark <- meta[-c(3,6),]
meta_noPark_RE_SMD <- metacont(Ne,Me,Se,Nc,Mc,Sc, data = meta_noPark,studlab =paste(Author),fixed = FALSE,random = TRUE,prediction = TRUE,sm="SMD",method.smd = "Cohen")
forest(meta_noPark_RE_SMD,digits.se = 1,digits.mean = 1,sortvar = TE,leftcols = c("Cancer type","Author"),leftlabs = c("Cancer type","Author"))

summary(meta_noPark_RE_SMD)


png(file = "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/FP_noPark.png", width = 3000, height = 2400, res = 300)
forest(meta_noPark_RE_SMD,digits.se = 1,digits.mean = 1,sortvar = TE,leftcols = c("Cancer type","Author"),leftlabs = c("Cancer type","Author"))
dev.off()


###### Funnel plot - no Park studies
png(file = "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/FunnelPlot_noLabels_NoPark.png", width = 2500, height = 2000, res = 300)
col.contour = c("lavender", "lightblue1", "lightcyan1")
funnel.meta(meta_noPark_RE_SMD, xlim = c(-0.6, 1),
            contour = c(0.9, 0.95, 0.99),
            col.contour = col.contour,studlab = FALSE)
legend(x = 0.6, y = 0.01, 
       legend = c("p < 0.1", "p < 0.05", "p < 0.01"),
       fill = col.contour)
dev.off()

#####################



######## Eggers test

eggers.test(meta_RE_SMD)
eggers.test(meta_noPark_RE_SMD)

metabias(meta_RE_SMD,method.bias = "linreg")

##################


################################# SENSITIVTY ANALYSIS - "leave one out analysis 
vars <- c("TE", "lower", "upper")
## setting up loop
metaOb <- metacont(Ne,Me,Se,Nc,Mc,Sc, data = meta[-c(1),],studlab =paste(Author),comb.fixed = TRUE,comb.random = TRUE,prediction = TRUE,sm="SMD",method.smd = "Cohen")
s1 <- summary(metaOb)
res.md.random <- data.frame(s1$random)[vars]
red.md.random <- round(res.md.random, 5)
names(res.md.random) <- c("Absolute difference", "CI lower", "CI upper")
summary_df_random <- res.md.random

for(i in 2:nrow(meta)){
  metaOb <- metacont(Ne,Me,Se,Nc,Mc,Sc, data = meta[-c(i),],studlab =paste(Author),comb.fixed = TRUE,comb.random = TRUE,prediction = TRUE,sm="SMD",method.smd = "Cohen")
  s1 <- summary(metaOb)
  res.md.random <- data.frame(s1$random)[vars]
  red.md.random <- round(red.md.random, 5)
  names(res.md.random) <- c("Absolute difference", "CI lower", "CI upper")
  summary_df_random <- rbind(summary_df_random,res.md.random)
}


row.names(summary_df_random) <- c(paste("Excluding:",meta$Author[seq(1:nrow(meta))]))
View(summary_df_random)


