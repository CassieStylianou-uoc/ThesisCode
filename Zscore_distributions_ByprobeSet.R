# Load the dataset from an RDS file
c <- readRDS("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2020_Directory/R_analysis/ProbeIterationScripts/ECAC_2021_QC/GWAS_DataFiles_2021/GWAS_2021_CNVSet_2_200.rds")

# Display the range of CNV lengths
range(c$length)

# Filter the dataset to show CNVs with length >= 3000000
c[c$length >= 3000000,]

# Create a contingency table of CNV calls (deletions and duplications) by cohort
table(c$cnv_call, c$cohort)

# Display the range of segment means for deletions
range(c$segmean[c$cnv_call == "del"]) # -8.99 to -3.00

# Display the range of segment means for duplications
range(c$segmean[c$cnv_call == "dup"]) # 2.00 to 5.98

# Display the range of BAF (B-allele frequency) scores, omitting missing values
range(na.omit(c$baf_score_negative))

# Plot the density of segment means for deletions with different minimum probe counts
plot(density(c$segmean[c$cnv_call == "del" & c$probes >= 2]), col = "red", ylim = c(0, 0.35), lwd = 2, main = "Z-score distribution across probe sets - deletions")
lines(density(c$segmean[c$cnv_call == "del" & c$probes >= 3]), col = "green", lwd = 2)
lines(density(c$segmean[c$cnv_call == "del" & c$probes >= 5]), col = "lightblue4", lwd = 2)
lines(density(c$segmean[c$cnv_call == "del" & c$probes >= 10]), col = "purple", lwd = 2)

# Plot the density of segment means for duplications with different minimum probe counts
plot(density(c$segmean[c$cnv_call == "dup" & c$probes >= 2]), col = "red", ylim = c(0, 0.7), lwd = 2, main = "Z-score distribution across probe sets - duplications")
lines(density(c$segmean[c$cnv_call == "dup" & c$probes >= 3]), col = "green", lwd = 2)
lines(density(c$segmean[c$cnv_call == "dup" & c$probes >= 5]), col = "lightblue4", lwd = 2)
lines(density(c$segmean[c$cnv_call == "dup" & c$probes >= 10]), col = "purple", lwd = 2)

# Extract Z-scores for deletions with different minimum probe counts
del2_Zscores <- c$segmean[c$probes >= 2 & c$cnv_call == "del"]
del3_Zscores <- c$segmean[c$probes >= 3 & c$cnv_call == "del"]
del5_Zscores <- c$segmean[c$probes >= 5 & c$cnv_call == "del"]
del10_Zscores <- c$segmean[c$probes >= 10 & c$cnv_call == "del"]

# Group the deletion Z-scores into a list
delZgroups <- list(del2_Zscores, del3_Zscores, del5_Zscores, del10_Zscores)

# Perform Kolmogorov-Smirnov (KS) tests for all pairwise comparisons of the deletion Z-score groups
for (i in 1:(length(delZgroups) - 1)) {
  for (j in (i + 1):length(delZgroups)) {
    ks_result <- ks.test(delZgroups[[i]], delZgroups[[j]])
    
    # Print the result of the KS test
    cat(sprintf("KS test between Group %d and Group %d:\n", i, j))
    print(ks_result)
    
    # Check the p-value to determine statistical significance
    if (ks_result$p.value < 0.05) {
      cat("The distributions are significantly different.\n\n")
    } else {
      cat("There is no significant difference in the distributions.\n\n")
    }
  }
}

# Extract Z-scores for duplications with different minimum probe counts
dup2_Zscores <- c$segmean[c$probes >= 2 & c$cnv_call == "dup"]
dup3_Zscores <- c$segmean[c$probes >= 3 & c$cnv_call == "dup"]
dup5_Zscores <- c$segmean[c$probes >= 5 & c$cnv_call == "dup"]
dup10_Zscores <- c$segmean[c$probes >= 10 & c$cnv_call == "dup"]

# Group the duplication Z-scores into a list
dupZgroups <- list(dup2_Zscores, dup3_Zscores, dup5_Zscores, dup10_Zscores)

# Perform Kolmogorov-Smirnov (KS) tests for all pairwise comparisons of the duplication Z-score groups
for (i in 1:(length(dupZgroups) - 1)) {
  for (j in (i + 1):length(dupZgroups)) {
    ks_result <- ks.test(dupZgroups[[i]], dupZgroups[[j]])
    
    # Print the result of the KS test
    cat(sprintf("KS test between Group %d and Group %d:\n", i, j))
    print(ks_result)
    
    # Check the p-value to determine statistical significance
    if (ks_result$p.value < 0.05) {
      cat("The distributions are significantly different.\n\n")
    } else {
      cat("There is no significant difference in the distributions.\n\n")
    }
  }
}