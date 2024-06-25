#### Reading in GWAS data (in xlsx format) 
delGWAS_3probe <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 4)
dupGWAS_3probe <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 5)
lofGWAS_3probe <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 6)

round(median(qchisq(delGWAS_3probe$p,df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.022
round(median(qchisq(dupGWAS_3probe$p,df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.022
round(median(qchisq(lofGWAS_3probe$p,df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.026

round(median(qchisq(delGWAS_3probe$p[delGWAS_3probe$All.overlaps>=6],df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.737
round(median(qchisq(dupGWAS_3probe$p[dupGWAS_3probe$All.overlaps>=6],df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.948
round(median(qchisq(lofGWAS_3probe$p[lofGWAS_3probe$All.overlaps>=6],df = 1,lower.tail = FALSE)) /qchisq(0.5,1),3) # 0.948




library(qqman)
alpha<-median(qchisq(1-delGWAS_3probe$p,1))/qchisq(0.5,1)
qq(delGWAS_3probe$p)
text(0.5,4, paste("lambda","=",  signif(alpha, digits = 3)) )


alpha<-median(qchisq(1-dupGWAS_3probe$p,1))/qchisq(0.5,1)
qq(dupGWAS_3probe$p)
text(0.5,4, paste("lambda","=",  signif(alpha, digits = 3)) )

alpha<-median(qchisq(1-lofGWAS_3probe$p,1))/qchisq(0.5,1)
qq(lofGWAS_3probe$p)
text(0.5,4, paste("lambda","=",  signif(alpha, digits = 3)))




##### Calculating lambda with only six 

alpha <- median(qchisq(1- delGWAS_3probe$p[delGWAS_3probe$All.overlaps>=6],1)/qchisq(0.5,1))
qq(delGWAS_3probe$p[delGWAS_3probe$All.overlaps>=6])
text(0.5,4, paste("lambda","=",  signif(alpha, digits = 3)) )


