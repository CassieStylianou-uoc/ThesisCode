
library(gdata)

PvalFunc <- function(df){
  ToreturnControl <- df %>% group_by(Control.overlaps) %>% summarise(min.p=min(p),max.p =max(p),min.case=min(Case.overlaps),max.case=max(Case.overlaps),x=n())
  ToreturnCase <- df %>% group_by(Case.overlaps) %>% summarise(min.p=min(p),max.p =max(p),min.control=min(Control.overlaps),max.control=max(Control.overlaps),x=n())
  ToreturnFull <- cbindX(ToreturnCase,ToreturnControl)
  return(ToreturnFull)
}


Pval_dels <- PvalFunc(del_threeProbe)
Pval_dups <- PvalFunc(dup_threeProbe)
Pval_lofs <- PvalFunc(lof_threeProbe)


PvalFuncAll <- function(df){
  Toreturn <- df %>% group_by(All.overlaps) %>% summarise(min.p=min(p),max.p=max(p),min.case=min(Case.overlaps),min.control=min(Control.overlaps),
                                                          max.case=max(Case.overlaps),max.control=max(Control.overlaps),x=n())
  return(Toreturn)
}


Pval_delsAll <- as.data.frame(PvalFuncAll(del_threeProbe))
Pval_dupsAll <- as.data.frame(PvalFuncAll(dup_threeProbe))
Pval_lofsAll <- as.data.frame(PvalFuncAll(lof_threeProbe))


library(xlsx)
write.xlsx(Pval_delsAll, file=paste0(dirP,"Pval_threeProbe.xlsx"), sheetName="DEL_Pval", row.names=FALSE)
write.xlsx(Pval_dupsAll, file=paste0(dirP,"Pval_threeProbe.xlsx"), sheetName="DUP_Pval", append=TRUE, row.names=FALSE)
write.xlsx(Pval_lofsAll, file=paste0(dirP,"Pval_threeProbe.xlsx"), sheetName="LOF_Pval", append=TRUE, row.names=FALSE)


