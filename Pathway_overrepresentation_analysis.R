######## Pathway enrichment 

library(gprofiler2)
library(readxl)

#### Reading in GWAS data (in xlsx format) 
delGWAS <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 4)
dupGWAS <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 5)
lofGWAS <- read_xlsx("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/GWAS_results/Compiled_GWASresults.xlsx",sheet = 6)

del_p <- delGWAS[delGWAS$p < 0.01,] #59
dup_p <- dupGWAS[dupGWAS$p < 0.01,] #58
lof_p <- lofGWAS[lofGWAS$p < 0.01,] #116 


dels_gost <- gprofiler2::gost(query = del_p$Gene.name,organism = "hsapiens",significant = T,evcodes = TRUE,sources = c("GO", "REAC", "HP", "HPA", "WP"),correction_method = "bonferroni" )
dels_results <- dels_gost$result

dups_gost <- gprofiler2::gost(query = dup_p$Gene.name,organism = "hsapiens",significant = T,evcodes = TRUE,sources = c("GO", "REAC", "HP", "HPA", "WP"),correction_method = "bonferroni" )
dups_results <- dups_gost$result


lofs_gost <- gprofiler2::gost(query = lof_p$Gene.name,organism = "hsapiens",significant = T,evcodes = TRUE,sources = c("GO", "REAC", "HP", "HPA", "WP"),correction_method = "bonferroni" )
lofs_results <- lofs_gost$result



########### Getting rid of any that are GO:MF/BP (MF=molecular function / BP = biological process and CC = cellular component)
####### https://biit.cs.ut.ee/gprofiler/page/docs

PA_delSS <- subset(dels_results,!grepl(c("GO:"),source)) #13
PA_dupSS <- subset(dups_results,!grepl("GO:",source)) # 1
PA_lofSS <- subset(lofs_results,!grepl("GO:",source)) #6

pathways <- append(PA_delSS$term_name,PA_dupSS$term_name)
pathways <- append(pathways,PA_lofSS$term_name) #20
unique_pathways <- unique(pathways) #15

pathDir <- "C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_3/ThreeProbeThesisResults/"
require(openxlsx)

list_of_datasets <- list("Deletion" = PA_delSS, "Duplication" = PA_dupSS,"Loss of function" = PA_lofSS,"Unique_pathways"=unique_pathways)
write.xlsx(list_of_datasets, file = paste0(pathDir,"PathwayCompiled.xlsx"))

