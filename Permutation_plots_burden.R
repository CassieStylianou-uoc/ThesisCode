


b <- read.csv("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/SampleSummary_ForBurden_3probes.csv")
b <- b [,-c(1)]


NumAll <- nrow(b)
NumCase <- nrow(b[b$case_control==1,])
NumControl <- nrow(b[b$case_control==0,])

PermFunction <- function(permDf,numCol,NumCase){
  i = 1
  b <- perm[,c(1,numCol)]
  x <- sample(b$OncID,NumCase,replace = T)
  b$Group <- ifelse(b$OncID %in% x,"G1","G2")
  M1 <- mean(b[,2][b$Group == "G1"])
  M2 <- mean(b[,2][b$Group == "G2"])
  t <- t.test(b[,2][b$Group=="G1"],b[,2][b$Group=="G2"],alternative = "two.sided",var.equal = FALSE)
  p <- t$p.value
  yy <- cbind(M1,M2,p,i)
  colnames(yy) <- c("G1.mean","G2.mean","p.value","Count") 
  df <- as.data.frame(yy)
  
  for(i in 2:10000){
    b <- perm[,c(1,numCol)]
    x <- sample(b$OncID,NumCase,replace = T)
    b$Group <- ifelse(b$OncID %in% x,"G1","G2")
    M1 <- mean(b[,2][b$Group == "G1"])
    M2 <- mean(b[,2][b$Group == "G2"])
    t <- t.test(b[,2][b$Group=="G1"],b[,2][b$Group=="G2"],alternative = "two.sided",var.equal = FALSE)
    p <- t$p.value
    yy <- cbind(M1,M2,p,i)
    colnames(yy) <- c("G1.mean","G2.mean","p.value","Count") 
    df <- rbind(df,yy)
  }
  for(i in 1:nrow(df)){
    df$logP[i] <- -log(df$p.value[i],2)
  }
  return(df)
}



perm_CNVs <- PermFunction(perm,2,NumCase)
perm_genicCNVs <- PermFunction(perm,3,NumCase)
perm_pLIcnvs <- PermFunction(perm,4,NumCase)
perm_WholepLICNVs <- PermFunction(perm,6,NumCase)
perm_ExonicCNVs <- PermFunction(perm,7,NumCase)
perm_CpGCNVs <- PermFunction(perm,8,NumCase)



makingPermPlots_onlyThree <- function(permutation_result,Three,y_limit){
  plot <- ggplot(data = permutation_result,aes(x=Count,y=logP))+geom_point(size = 0.5)+theme_classic()+
    geom_hline(yintercept = -log(Three,2),color="black",linewidth = 1)+
    scale_x_continuous(expand = c(0,0),limits=c(0,10050)) + scale_y_continuous(expand = c(0,0),limits = c(0,y_limit)) + xlab(label = "Permutation count") +
    ylab(label = "-log10(P-value)") 
}



library(ggplot2)
p_cnvPlot <- makingPermPlots_onlyThree(perm_CNVs,4.37E-63,240)
p_genePlot <- makingPermPlots_onlyThree(perm_genicCNVs,2.10E-50,200)
p_pliPlot <- makingPermPlots_onlyThree(perm_pLIcnvs,3.58E-21,80)
p_WholepLIplot <- makingPermPlots_onlyThree(perm_WholepLICNVs,0.000182,25)
p_exonicCNVs <- makingPermPlots_onlyThree(perm_ExonicCNVs,2.43E-41,180)
p_CpGPlot <- makingPermPlots_onlyThree(perm_CpGCNVs,1.16E-27,100)




dir <- c("C:/Users/styca196/OneDrive - University of Otago/cstylianou/2023/Results_2/final_ECAC_burden/")




jpeg(paste0(dir,"permCNVPlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_cnvPlot
dev.off()

jpeg(paste0(dir,"permGenicPlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_genePlot
dev.off()

jpeg(paste0(dir,"permPLIPlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_pliPlot
dev.off()

jpeg(paste0(dir,"permPLIWholePlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_WholepLIplot
dev.off()

jpeg(paste0(dir,"permExonicPlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_exonicCNVs
dev.off()

jpeg(paste0(dir,"permCpGPlot_onlyThree_replacement.jpg"), width = 240, height = 300)
p_CpGPlot
dev.off()

