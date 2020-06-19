  qq=DEG[DEG$change !='NOT',]
nrow(DEG[DEG$change !='NOT',])

diffgene=DEG[DEG$change !='NOT',]
write.csv(diffgene,"diffgene.csv")
diffgene_exprset=exprSet[rownames(diffgene),]
two_in_one_table=cbind(diffgene,diffgene_exprset)
write.csv(diffgene_exprset,"diffgene_exprset.csv")
write.csv(two_in_one_table,"two_in_one_table.csv")







DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
                       qq=DEG[DEG$change !='NOT',]                       
                       
write.csv(DEG,"entrezid_dif_Calcified_Normal.csv")
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')