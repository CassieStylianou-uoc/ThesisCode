
###### Binary matrixes made in MakingBinaryMatrixes.R file 



x = list.files(pattern = "Probe_")
outfile = gsub(".rds","_association.rds",x)
Samples = readRDS("GWAS_2021_SampleSet.rds")


for(j in 1:length(x)){
  tmp <- readRDS(x[j]) ## binary matrix
  Samples_tmp <- rownames(tmp) #21,933 samples 
  tmp_SamplesCohort <- Samples[Samples$OncID %in% Samples_tmp,]
  tmp_SamplesCohort <- tmp_SamplesCohort[match(rownames(tmp),tmp_SamplesCohort$OncID),]
  
  case.contol <- tmp_SamplesCohort$case_control
  ass_tmp=apply(tmp,2, function(i) 
    glm(tmp_SamplesCohort$case_control ~ i, family=binomial))
  r=lapply(ass_tmp, function(i) 
    data.frame(OR=exp(i$coefficients[2]), p=coef(summary(i))[2,4], CI=paste(exp(confint(i))[2,],collapse ="-"))
  )
  final.table = data.table::rbindlist(r, idcol = "name")
  
  saveRDS(final.table,file=outfile[j])
}


