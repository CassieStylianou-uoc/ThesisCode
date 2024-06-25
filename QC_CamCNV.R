### Reading in original sample level data obtained from Joe Dennis 
og_samples <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/CNV_data_files/ECAC_BCAC_CNVData/Original_data/ecac_samples.csv")
og_clin <- read.delim("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/CNV_data_files/ECAC_BCAC_CNVData/Original_data/ecac_oncoarray_logitregress_may2016.txt")

# DLRS exclusiso: Filtering samples based on DLRS values before and after PCA
samples_DLRS_pre <- og_samples[og_samples$DLRS_before_PCA <=0.2,]
table(samples_DLRS_pre$case_control) #18,242 / 4285
samples_DLRS_pre[is.na(samples_DLRS_pre)] <-0
samples_DLRS_post <- samples_DLRS_pre[samples_DLRS_pre$DLRS_after_PCA <=0.1,]
table(samples_DLRS_post$case_control) #17.820 / 4.218
DLRS_passingSamples <- samples_DLRS_post$OncID

## Read in orginal CNV data - split by CNV type 
og_dels<- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/CNV_data_files/ECAC_BCAC_CNVData/Original_data/ecac_dels_pre_qc.csv")
og_dups<- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/CNV_data_files/ECAC_BCAC_CNVData/Original_data/ecac_dups_pre_qc.csv")

# Commbine CNVs into a single dataframe 
og_cnvs <- rbind(og_dels,og_dups)

# Annotate CNV data with case/control cohort information
for(i in 1:nrow(og_cnvs)){
  o <- og_samples[og_samples$OncID == og_cnvs$sampleid[i],] ## Match CNV sample ID with sample data
  og_cnvs$cohort[i] <- o$case_control # Assign case/control status
  print(i)
}

### Pre-QC: Display initial CNV call distribution by cohort
table(og_cnvs$cohort,og_cnvs$cnv_call)

### Filter CNVs based on samples that passed DLRS exclusions
cnvs_passingDLRS <- og_cnvs[og_cnvs$sampleid %in% DLRS_passingSamples,]
####### # Display CNV QC exclusions
table(cnvs_passingDLRS$cnv_call[cnvs_passingDLRS$qc_category == ""],cnvs_passingDLRS$cohort[cnvs_passingDLRS$qc_category == ""])

# CNV exclusions by cohort: Display CNV calls by QC category for cases and controls
table(cnvs_passingDLRS$cnv_call[cnvs_passingDLRS$cohort == 1],cnvs_passingDLRS$qc_category[cnvs_passingDLRS$cohort == 1])
table(cnvs_passingDLRS$cnv_call[cnvs_passingDLRS$cohort == 0],cnvs_passingDLRS$qc_category[cnvs_passingDLRS$cohort == 0])

# Annotate CNVs with study information
for(i in 1:nrow(cnvs_passingDLRS)){
  x <- og_samples$ecac_study[og_samples$OncID == cnvs_passingDLRS$sampleid[i]]
  cnvs_passingDLRS$ecac_study[i] <- x
  print(i)
}
table(cnvs_passingDLRS$ecac_study,cnvs_passingDLRS$qc_category)

# Probe range exclusions: Display CNV call distribution for CNVs with <= 200 probes by cohort
table(cnvs_passingDLRS$cnv_call[cnvs_passingDLRS$qc_category == "" & cnvs_passingDLRS$probes <=200],cnvs_passingDLRS$cohort[cnvs_passingDLRS$qc_category == "" & cnvs_passingDLRS$probes <=200])

### CNV load exclusions: Tally CNV counts per sample based on minimum probe counts
sample_cnv_tally <- as.data.frame(DLRS_passingSamples)
colnames(sample_cnv_tally) <- c("OncID")
for(i in 1:nrow(sample_cnv_tally)){
  cnvs <- cnvs_passingDLRS_passingJoeQC_probeRange[cnvs_passingDLRS_passingJoeQC_probeRange$sampleid == sample_cnv_tally$OncID[i],]
  sample_cnv_tally$cnv_2min[i] <- nrow(cnvs)
  sample_cnv_tally$cnv_3min[i] <- nrow(cnvs[cnvs$probes >=3,])
  sample_cnv_tally$cnv_4min[i] <- nrow(cnvs[cnvs$probes >=4,])
  sample_cnv_tally$cnv_5min[i] <- nrow(cnvs[cnvs$probes >=5,])
  sample_cnv_tally$cnv_6min[i] <- nrow(cnvs[cnvs$probes >=6,])
  sample_cnv_tally$cnv_7min[i] <- nrow(cnvs[cnvs$probes >=7,])
  sample_cnv_tally$cnv_8min[i] <- nrow(cnvs[cnvs$probes >=8,])
  sample_cnv_tally$cnv_9min[i] <- nrow(cnvs[cnvs$probes >=9,])
  sample_cnv_tally$cnv_10min[i] <- nrow(cnvs[cnvs$probes >=10,])
  print(i)
}

# Annotate samples with clinical and demographic information
postQC_samples_globalCNV <- sample_cnv_tally
for(i in 1:nrow(postQC_samples_globalCNV)){
  info <- og_samples[og_samples$OncID == postQC_samples_globalCNV$OncID[i],]
  clin_info <- og_clin[og_clin$id == postQC_samples_globalCNV$OncID[i],]
  postQC_samples_globalCNV$ecac_study[i] <- info$ecac_study
  postQC_samples_globalCNV$ecac_strata[i] <- info$ecac_strata
  postQC_samples_globalCNV$case_control[i] <- info$case_control
  postQC_samples_globalCNV$my_status_clin[i] <- clin_info$mystatus
  postQC_samples_globalCNV$hist[i] <- clin_info$hist
  postQC_samples_globalCNV$age[i] <- clin_info$age
  postQC_samples_globalCNV$weight[i] <- clin_info$weight
  postQC_samples_globalCNV$bmi[i] <- clin_info$bmi
  print(i)
}
# Exclude samples with more than 50 CNVs
max50_CNV <- postQC_samples_globalCNV[postQC_samples_globalCNV$cnv_2min <=50,]
table(max50_CNV$case_control) #17,818 / 4,115
