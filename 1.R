rm(list=ls())
library("affy")
library("limma")
#library("affyPLM")
library("RColorBrewer")
library("pheatmap")
library("affycoretools")
#library("simpleaffy")
library("BiocGenerics")
library("GEOquery")
#53986
eSet <- getGEO('GSE119717', destdir=".",AnnotGPL = F,getGPL = F)
save(eSet,file='1.1eSet.Rdata')

rm(list=ls())
load('1.1eSet.Rdata')
b = eSet[[1]]
raw_exprSet=exprs(b) 
raw_exprSet[1:4,1:4]
# > raw_exprSet[1:4,1:4]
# GSM1304836 GSM1304837 GSM1304838 GSM1304839
# 1415670_at   1574.460   1432.630   1328.000   1171.610
# 1415671_at  10303.400  10683.800  10038.300  10209.600
# 1415672_at   9121.120   8679.330   9367.220   9710.290
# 1415673_at    596.757    501.811    420.459    445.928
exprSet=log2(raw_exprSet)
exprSet[1:4,1:4]
# > exprSet[1:4,1:4]
# GSM1304836 GSM1304837 GSM1304838 GSM1304839
# 1415670_at   10.62064   10.48445  10.375039  10.194277
# 1415671_at   13.33083   13.38314  13.293227  13.317639
# 1415672_at   13.15500   13.08337  13.193405  13.245299
# 1415673_at    9.22100    8.97100   8.715821   8.800667

# IF   IF   IF 
#exprSet=raw_exprSet

phe=pData(b)
library(stringr)
group_list= str_split(as.character(phe$title),',',simplify = T)[,2]
#group_list= str_split(as.character(phe$title),'_',simplify = T)[,1]
#group_list= as.character(phe[,34])
#= as.character(phe$`high surfactant (yes or no):ch1`)
# group_list= str_split(as.character(phe$title),' ',simplify = T)[,3]
# save(raw_exprSet,group_list,phe,ids,
#      file='GSE8021_raw_exprSet_IDO.Rdata')
# phe$title
# group_list= str_split(as.character(phe$title),' ',simplify = T)[,11]
# group_list=paste0('group',group_list,'h')
# group_list=as.character(datTraits$ratio)

#check  the  type 
#> group_list
#[1] " untreated" " untreated" " untreated" " untreated" " IFNg"      " IFNg"      " IFNg"      " IFNg"     
#[9] " LPS"       " LPS"       " LPS"       " LPS"       " IFNg+LPS"  " IFNg+LPS"  " IFNg+LPS"  " IFNg+LPS" 
save(exprSet,group_list,file='1.2exprSet_group_list.Rdata')

rm(list=ls())
load('1.2exprSet_group_list.Rdata')
# library("hgu133a2.db")
# ids=toTable(mouse4302SYMBOL)
# ids=toTable(hgu133a2SYMBOL)
# ids <- read.delim("GPL5356.csv")
# ids=ids[,c(1,3)]
ids <- read.csv("~/Desktop/XIUGAI/GPL3.csv")
colnames(ids)=c("probe_id","symbol")  
head(ids)
length(unique(ids$symbol))
tail(sort(table(ids$symbol)))
table(sort(table(ids$symbol)))
plot(table(sort(table(ids$symbol))))
#判断前面一个向量内的元素是否在后面一个向量中，返回布尔值
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]

jimmy <- function(exprSet,ids){
  tmp = by(exprSet,
           ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}
new_exprSet <- jimmy(exprSet,ids)
save(new_exprSet,ids,group_list,file='1.3gene_symbol_exprSet_group_list.Rdata')

#load('1.2exprSet_group_list.Rdata')
rm(list=ls())
load(file='1.3gene_symbol_exprSet_group_list.Rdata')
# new_exprSet=diffgene_exprset[hub$x,]
# colnames(new_exprSet)=group_list
# write.csv(new_exprSet,"93sample_93hubgene.csv")
# exprSet=new_exprSet
#group_list=as.character(datTraits[,2])

if(T){dat=new_exprSet#每次都要检测数据
dat=t(dat)#画PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
dat=as.data.frame(dat)#将matrix转换为data.frame
dat=cbind(dat,group_list) #cbind横向追加，即将分组信息追加到最后一列
library("FactoMineR")#画主成分分析图需要加载这两个包
library("factoextra") # The variable group_list (index = 54676) is removed # before PCA analysis
#现在dat最后一列是group_list，需要重新赋值给一个dat.pca,这个矩阵是不含有分组信息的
dat.pca <- PCA(dat[,-ncol(dat)], graph = FALSE)
fviz_pca_ind(dat.pca,geom.ind = "point", # show points only (nbut not "text")
             col.ind = dat$group_list, # color by groups # palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, legend.title = "Groups")# Concentration ellipses
#ggsave('2.1_PCA.png')
}

if(T){  ## hclust
  exprSet=new_exprSet
  colnames(exprSet)=paste(group_list,1:ncol(exprSet),sep='_')
  # Define nodePar
  nodePar <- list(lab.cex = 0.6, pch = c(NA, 19),cex = 0.7, col = "blue")
  hc=hclust(dist(t(exprSet)))
  par(mar=c(5,5,5,10))
  plot(as.dendrogram(hc), nodePar = nodePar, horiz = TRUE)
  ## PCA
  library(ggfortify)
  df=as.data.frame(t(exprSet))
  df$group=group_list
  autoplot(prcomp(df[,1:(ncol(df)-1)]), data=df,colour = 'group') + theme_bw()}

if(T){
  library(reshape2)
  exprSet_L=melt(exprSet)
  colnames(exprSet_L)=c('probe','sample','value')
  exprSet_L$group=rep(group_list,each=nrow(exprSet))
  head(exprSet_L)
  library(ggplot2)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_violin()
  print(p)
  p=ggplot(exprSet_L,aes(value,fill=group))+geom_histogram(bins = 200)+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()+facet_wrap(~sample, nrow = 4)
  print(p)
  p=ggplot(exprSet_L,aes(value,col=group))+geom_density()
  print(p)
  p=ggplot(exprSet_L,aes(x=sample,y=value,fill=group))+geom_boxplot()
  p=p+stat_summary(fun.y="mean",geom="point",shape=23,size=3,fill="red")
  p=p+theme_set(theme_set(theme_bw(base_size=20)))
  p=p+theme(text=element_text(face='bold'),axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank())
  print(p) 
}
