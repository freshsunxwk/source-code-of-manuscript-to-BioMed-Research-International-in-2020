rm(list=ls())
load(file='1.3gene_symbol_exprSet_group_list.Rdata')
exprSet=new_exprSet
#colnames(exprSet)
#exprSet=NOcsv
#exprSet=exprSet[,-1]
dim(exprSet)
group_list
#group_list=c("LPS","LPS","LPS","WT","WT","WT","LPS","LPS","LPS","WT","WT","WT")
# exprSet=exprSet[,c(1:4,9:12)]
# group_list=group_list[c(1:4,9:12)]
# group_list=c("WT","WT","WT","WT","LPS","LPS","LPS","LPS")                         
library(limma)
# tmp=data.frame(case=c(0,0,0,1,1,1),control=c(1,1,1,0,0,0))
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design
# contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
# contrast.matrix<-makeContrasts("LPS-WT",levels = design)
# contrast.matrix<-makeContrasts("TAA-Contro",levels = design)
contrast.matrix<-makeContrasts("yes-no",levels = design)
contrast.matrix<-makeContrasts("DPVAT-NDPVAT",levels = design)
DPVAT NDPVAT
contrast.matrix ##这个矩阵声明，我们要把72h处理组跟未处理的细胞系进行差异分析比较
save(exprSet,design,contrast.matrix,file="1.4forDEG.Rdata")

rm(list=ls())
load(file='1.4forDEG.Rdata')
deg = function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)##step1#
  fit2 <- contrasts.fit(fit, contrast.matrix) #step2##这一步很重要，大家可以自行看看效果
  fit2 <- eBayes(fit2)  ## default no trend !!!##eBayes() with trend=TRUE #
  tempOutput = topTable(fit2, coef=1, n=Inf)#step3
  nrDEG = na.omit(tempOutput) 
  #write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
  head(nrDEG)
  return(nrDEG)}
re = deg(exprSet,design,contrast.matrix)
nrDEG=re

library(pheatmap)## heatmap
choose_gene=head(rownames(nrDEG),50) ##50 maybe better
choose_matrix=exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,filename = 'DEG_top50_heatmap.png')

library(ggplot2)## volcano plot

colnames(nrDEG)
plot(nrDEG$logFC,-log10(nrDEG$P.Value))
DEG=nrDEG
#logFC_cutoff <- with(DEG,mean(abs( logFC)) + 2*sd(abs( logFC)))
# logFC_cutoff=1
#logFC_cutoff=0.585
logFC_cutoff=0.25
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',]))
g = ggplot(data=DEG,aes(x=logFC, y=-log10(P.Value),color=change))+geom_point(alpha=0.5, size=1.20) +
  theme_set(theme_set(theme_bw(base_size=20)))+xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+xlim(-2,2)+
  scale_colour_manual(values = c('blue','black','red')) ## corresponding to the levels(res$change)
print(g)
ggsave(g,filename = 'volcano.png')
write.csv(DEG,"DEG.csv")
save(exprSet,group_list,nrDEG,DEG, file='1.5after_DEG.Rdata')
  
rm(list=ls())
load(file='1.5after_DEG.Rdata')
library(clusterProfiler)
source('functions.R')
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(rownames(DEG), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
head(DEG)
DEG$SYMBOL = rownames(DEG)
DEG=merge(DEG,df,by='SYMBOL')
head(DEG)

gene_up= DEG[DEG$change == 'UP','ENTREZID'] 
gene_down=DEG[DEG$change == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )

data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

## KEGG pathway analysis
if(T){
  ###   over-representation test
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.05)
  head(kk.down)[,1:6]
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  
  
  write.csv(kegg_diff_dt,"kegg_diff_dt.csv")
  write.csv(kegg_down_dt,"kegg_down_dt.csv")
  write.csv(kegg_up_dt,"kegg_up_dt.csv")
  
  
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  
  ggsave(g_kegg,filename = 'kegg_up_down.png')
  
  ###  GSEA 
  kk_gse <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 120,
                    pvalueCutoff = 0.9,
                    verbose      = FALSE)
  head(kk_gse)[,1:6]
  gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
  
  down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
  up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
  
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down_gsea.png')
  
  #write.csv(kegg_up_dt,"kegg_up_dt.csv")
}

### GO database analysis 

g_list=list(gene_up=gene_up,
            gene_down=gene_down,
            gene_diff=gene_diff)

if(F){
  go_enrich_results <- lapply( g_list , function(gene) {
    lapply( c('BP','MF','CC') , function(ont) {
      cat(paste('Now process ',ont ))
      ego <- enrichGO(gene          = gene,
                      universe      = gene_all,
                      OrgDb         = org.Hs.eg.db,
                      ont           = ont ,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.99,
                      qvalueCutoff  = 0.99,
                      readable      = TRUE)
      
      print( head(ego) )
      return(ego)
    })
  })
  save(go_enrich_results,file = 'go_enrich_results.Rdata')
  
}


load(file = 'go_enrich_results.Rdata')

n1= c('gene_up','gene_down','gene_diff')
n2= c('BP','MF','CC') 
for (i in 1:3){
  for (j in 1:3){
    fn=paste0('dotplot_',n1[i],'_',n2[j],'.png')
    cat(paste0(fn,'\n'))
    png(fn,res=150,width = 1080)
    print( dotplot(go_enrich_results[[i]][[j]] ))
    dev.off()
  }
}

