datTraits=datTraitst
  datTraits = data.frame(GSM=phe[,2],
                         outcome=group_list,ratio=phe[,33])
  datTraits = data.frame(GSM=phe[,2],
                         ratio=group_list)
  
  rownames(datTraits)=datTraits$GSM
  write.csv(datTraits,"datTraitsto010101.csv")
  save(datTraits,file='datTraits.Rdata')
  
  
  datTraits <- read.csv("D:/DEG/1129again/datTraits-degree.csv", row.names=1)

  #fpkm=diffgene_exprset
  fpkm=new_exprSet
  
  
  library(GEOquery)
  a=getGEO('GSE48213')
  metadata=pData(a[[1]])[,c(2,10,12)]
  datTraits = data.frame(gsm=metadata[,1],
             cellline=trimws(sapply(as.character(metadata$characteristics_ch1),function(x) strsplit(x,":")[[1]][2])),
             subtype=trimws(sapply(as.character(metadata$characteristics_ch1.2),function(x) strsplit(x,":")[[1]][2]))
             )
  save(fpkm,datTraits,file = 'GSE9106-wgcna-input-new8.RData')
}
load('GSE52093-wgcna-input.RData')
library(WGCNA)
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=4

## step 1 :
datTraits=datTraits[-c(67,87,83),] 
fpkm =fpkm[,-c(67,87,83)] 
datTraits[,2]=="testing set"
  fpkm[1:4,1:4]
  head(datTraits)
  dim(fpkm)
  m.mad <- apply(new_exprSet,1,mad)
  new_exprSet_cut <- new_exprSet[which(m.mad > 
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
  RNAseq_voom=new_exprSet_cut
  
  fpkm=10^(fpkm)
  #boxplot(fpkm)
 
  exprSet[1:4,1:4]
  exprSet
  
  RNAseq_voom <- fpkm 
  ## 因为WGCNA针对的是基因进行聚类，而一般我们的聚类是针对样本用hclust即可，所以这个时候需要转置
  #WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
  WGCNA_matrix = t(RNAseq_voom[order(apply(RNAseq_voom,1,mad), decreasing = T)[1:5000],])
  datExpr0 <- WGCNA_matrix  ## top 5000 mad genes
  datExpr <- datExpr0 
  
  ## 下面主要是为了防止临床表型与样本名字对不上
  
  #rownames(datExpr) <- gsub(".CEL.gz", "", rownames(datExpr))
  
  sampleNames = rownames(datExpr);
  traitRows = match(sampleNames, datTraits$GSM)
  rownames(datTraits) = datTraits[traitRows, 1]
  

  #pdf(file="Fig.pdf",width=80, height=60, units='mm',res=1200,compression='lzw',pointsize=6)

  ## step 2 
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  pdf("step2-beta-value.pdf")
  # Plot the results:
  ##sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()

结果
sft$powerEstimate=14
12


pdf("step2-Scale-free.pdf")
k <- softConnectivity(datE=datExpr,power=sft$powerEstimate) 
sizeGrWindow(10, 5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()
## step3 构建加权共表达网络（Weight co-expression network)
## 首先是一步法完成网络构建

    net = blockwiseModules(
      datExpr,
      power = sft$powerEstimate,
      maxBlockSize = 6000,
      TOMType = "unsigned", minModuleSize = 30,
      reassignThreshold = 0, mergeCutHeight = 0.25,
      numericLabels = TRUE, pamRespectsDendro = FALSE,
      saveTOMs = TRUE,
      saveTOMFileBase = "save-TOM",
      verbose = 3
    )
    table(net$colors) 
 # 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14 
 # 2735 2635 1727  712  650  517  252  219  198  111   63   57   46   42   36
  
sft$powerEstimate =6
折叠起来
## 然后是分布法完成网络构建，仅供有探索精神的同学挑战。######

## 构建加权共表达网络分为两步：
## 1. 计算邻近值，也是就是两个基因在不样品中表达量的表达相关系数(pearson correlation rho)，
## 参考 2.b.2 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
## 2. 计算topology overlap similarity (TOM)。 WGCNA认为，只通过计算两个基因的表达相关系数构建共表达网络是不足够的。
## 于是他们用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
## 参考 2.b.3 in https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf


save(TOM,file = 'TOM.RData')#

	#(1)网络构建 Co-expression similarity and adjacency 
	adjacency = adjacency(datExpr, power = sft$powerEstimate) 
	#(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
	TOM = TOMsimilarity(adjacency);
	
	dissTOM = 1-TOM
	# (3) 聚类拓扑矩阵 Call the hierarchical clustering function
	geneTree = hclust(as.dist(dissTOM), method = "average");
	
	
	# Plot the resulting clustering tree (dendrogram)
	sizeGrWindow(12,9)
	## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
	plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
		 labels = FALSE, hang = 0.04);
	#(4) 聚类分支的修整 dynamicTreeCut 
	# We like large modules, so we set the minimum module size relatively high:
	minModuleSize = 30;
	# Module identification using dynamic tree cut:
	dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
								deepSplit = 2, pamRespectsDendro = FALSE,
								minClusterSize = minModuleSize);
	table(dynamicMods)
	#dynamicMods
#	1    2    3    4    5    6    7    8 
#	4836 4204  379  168  148   95   88   82 
	#4. 绘画结果展示
	# Convert numeric lables into colors
	dynamicColors = labels2colors(dynamicMods)
	table(dynamicColors)
#	dynamicColors
#	black      blue     brown     green      pink       red turquoise    yellow 
#	88      4204       379       148        82        95      4836       168 
	# Plot the dendrogram and colors underneath
	#sizeGrWindow(8,6)
	plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
						dendroLabels = FALSE, hang = 0.03,
						addGuide = TRUE, guideHang = 0.05,
						main = "Gene dendrogram and module colors")
	# (5)聚类结果相似模块的融合，Merging of modules whose expression profiles are very similar
	#在聚类树中每一leaf是一个短线，代表一个基因，
	#不同分之间靠的越近表示有高的共表达基因，将共表达极其相似的modules进行融合
	# Calculate eigengenes
	MEList = moduleEigengenes(datExpr, colors = dynamicColors)
	MEs = MEList$eigengenes
	# Calculate dissimilarity of module eigengenes
	MEDiss = 1-cor(MEs);
	# Cluster module eigengenes
	METree = hclust(as.dist(MEDiss), method = "average");
	# Plot the result
	
	
	
	#pdf("Clustering-of-module-eigengenes.pdf",width=50, height=60, units='mm',res=1200,compression='lzw')
	sizeGrWindow(7, 6)
	plot(METree, main = "Clustering of module eigengenes",
		 xlab = "", sub = "")
	#选择有75%相关性的进行融合
	MEDissThres = 0.25
	# Plot the cut line into the dendrogram
	abline(h=MEDissThres, col = "red")
#	dev.off()
	
	
	# Call an automatic merging function
	merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
	# The merged module colors
	mergedColors = merge$colors;
	# Eigengenes of the new merged modules:
	mergedMEs = merge$newMEs
	
	table(mergedColors)
#	mergedColors
#	black         blue        brown         cyan        green  greenyellow         grey 
	255         1403          361           83          326          201            6 
	grey60   lightgreen  lightyellow midnightblue       purple          red       salmon 
	61           45           42           77          483          259           84 
	tan    turquoise 
	101         1213 
	
	pdf("Dynamic-Tree-Cut-Merged-dynamic.pdf",width=90, height=60, units='mm',res=1200,compression='lzw',pointsize=6)
	
	pdf("Dynamic-Tree-Cut-Merged-dynamic.pdf")
	plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
	                    c("Dynamic Tree Cut", "Merged dynamic"),
	                    dendroLabels = FALSE, hang = 0.03,
	                    addGuide = TRUE, guideHang = 0.05)
	dev.off()

#####
	plotEigengeneNetworks(mergedMEs, 
	                      "Eigengene adjacency heatmap", 
	                      marHeatmap = c(3,4,2,2), 
	                      plotDendrograms = FALSE, 
	                      xLabelsAngle = 90)
	###############################
	# 画出指定模块表达量的热图：     not  work  lllll
	which.module=Freq_MS_max$GS_color; 
	ME=mergedMEs[, paste("ME",which.module, sep="")]
	par(mfrow=c(2,1), mar=c(0,4.1,4,2.05))
	plotMat(t(scale(multiExpr[[1]]$data[,colorh1==which.module ]) ),
	        nrgcols=30,rlabels=F,rcols=which.module,
	        main=which.module, cex.main=2)
	par(mar=c(2,2.3,0.5,0.8))
	barplot(ME, col=which.module, main="", cex.main=2,
	        ylab="eigengene expression",xlab="array sample")
	
	
	datTraitsbeiyong=datTraits
	# 画出样本聚类图（上）与样本性状热图（下）：
	str(datTraits)
	datTraits=datTraits[,-1]
	traitColors = numbers2colors(datTraits, signed = TRUE,centered=TRUE);
	datExpr_tree<-hclust(dist(datExpr), method = "average")
	pdf("sample-subtype-cluster.pdf")
	plotDendroAndColors(datExpr_tree, 
	                    traitColors, 
	                    groupLabels = names(datTraits), 
	                    rowTextAlignment = "right-justified",
	                    addTextGuide = TRUE ,
	                    hang = 0.03,
	                    dendroLabels = NULL, # 是否显示树labels
	                    addGuide = FALSE,  # 显示虚线
	                    guideHang = 0.05,
	                    cex.dendroLabels = 0.5,
	                    cex.rowText = 0.5,
	                    main = "Sample dendrogram and trait heatmap")
	dev.off()
	nGenes = ncol(datExpr)
	nSamples = nrow(datExpr)
	
	moduleTraitCor_noFP <- cor(mergedMEs, datTraits[,1:5], use = "p");
	moduleTraitPvalue_noFP = corPvalueStudent(moduleTraitCor_noFP, nSamples); 
	textMatrix_noFP <- paste(signif(moduleTraitCor_noFP, 2), "\n(", signif(moduleTraitPvalue_noFP, 1), ")", sep = ""); 
	pdf("step5-Module-trait-relationships.pdf")
	
	par(mar = c(10, 8.5, 3, 3)); 
	labeledHeatmap(Matrix = moduleTraitCor_noFP, 
	               xLabels = names(datTraits[,1:5]), 
	               yLabels = names(mergedMEs), 
	               ySymbols = names(mergedMEs), 
	               colorLabels = FALSE, 
	               colors = blueWhiteRed(50), 
	               textMatrix = textMatrix_noFP,
	               setStdMargins = FALSE, 
	               cex.text = 0.65, 
	               zlim = c(-1,1), 
	               main = paste("Module-trait relationships"))
	dev.off()
	
	
	nGenes = ncol(datExpr)
	nSamples = nrow(datExpr)
	cor_ADR <- signif(WGCNA::cor(datTraits,mergedMEs,use="p",method="pearson"),5)
	p.values <- corPvalueStudent(cor_ADR,nSamples=nrow(datTraits))
	Freq_MS_max_cor <- which.max(abs(cor_ADR))
	Freq_MS_max_p <- which.min(p.values)
	
	
	
	GS1 <- as.numeric(WGCNA::cor(datTraits[,1],datExpr,use="p",method="pearson"))
	# 显著性是绝对值：
	GeneSignificance <- abs(GS1)
	# 获得该性状在每个模块中的显著性：
	ModuleSignificance <- tapply(GeneSignificance,mergedColors,mean,na.rm=T)
which.max(ModuleSignificance)
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,mergedColors,ylim=c(-0,0.6))

#Barplot of module significance defined as the mean gene significance across all genes in the module. The
#green and brown modules are the most promising.



softPower=6
# 计算每个基因模块内部连接度，也就是基因直接两两加权相关性。
ADJ1=abs(cor(datExpr,use="p"))^softPower 
# 根据上面结果和基因所属模块信息获得连接度：
# 整体连接度 kTotal，模块内部连接度：kWithin，kOut=kTotal-kWithin， kDiff=kIn-kOut=2*kIN-kTotal 
#colorh1=turquoise#
Alldegrees1=intramodularConnectivity(ADJ1, mergedColors) 

# 注意模块内基于特征向量基因连接度评估模块内其他基因： de ne a module eigengene-based connectivity measure for each gene as the correlation between a the gene expression and the module eigengene
# 如 brown 模块内：kM Ebrown(i) = cor(xi, MEbrown) ， xi is the gene expression pro le of gene i and M Ebrown is the module eigengene of the brown module
# 而 module membership 与内部连接度不同。MM 衡量了基因在全局网络中的位置。
datKME=signedKME(multiExpr[[1]]$data, datME, outputColumnName="MM.")#
######################################mei pao tong 
datKME  #              not work kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk

mergedColors
black         blue        brown         cyan        green  greenyellow         grey 
255         1403          361           83          326          201            6 
grey60   lightgreen  lightyellow midnightblue       purple          red       salmon 
61           45           42           77          483          259           84 
tan    turquoise 
101         1213 


which.color=Freq_MS_max$GS_color; #
which.color="grey60"
restrictGenes=colorh1==which.color #
restrictGenes=mergedColors==which.color 
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes], 
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^4, 
                   col=which.color, 
                   xlab="Intramodular Connectivity", 
                   ylab="(Module Membership)^4")



GS_spe <-  as.numeric(apply(datExpr,2,function(x){
  biserial.cor(x,datTraits[,1])
}))

GS_spe <-  as.numeric(apply(datExpr,2,function(x){
  biserial.cor(x,datTraits[,trait_minP$trait[i]])
}))
GeneSignificance_spe <- abs(GS_spe)
# 基于显著性和MM计算每个基因与 指定trait 的关联，结果包括p, q, cor, z, 
NS1=networkScreening(y=allTraits[,trait_minP$trait[[i]]], 
                     datME=datME, 
                     datExpr=multiExpr[[1]]$data, 
                     oddPower=3, 
                     blockSize=1000, 
                     minimumSampleSize=4, 
                     addMEy=TRUE, 
                     removeDiag=FALSE, 
                     weightESy=0.5) 

NS1=networkScreening(y=datTraits[,1], 
                     datME=MEs, 
                     datExpr=datExpr, 
                     oddPower=3, 
                     blockSize=1000, 
                     minimumSampleSize=4, 
                     addMEy=TRUE, 
                     removeDiag=FALSE, 
                     weightESy=0.5) 
rownames(NS1) <- colnames(datExpr)
kME_MM.turquoise
FilterGenes_spe = ((GeneSignificance> 0.2) & (abs(datKME$kME_MM.turquoise)>0.8) & (NS1$q.Weighted < 0.01) ) 

# 根据 基因与指定性状的直接相关性(biserial.cor)，模块身份，和加权相关性 筛选基因：
FilterGenes_spe = ((GeneSignificance_spe > 0.2) & (abs(datKME[paste("MM.",Freq_MS_max$GS_color,sep="")])>0.8) & (NS1$q.Weighted < 0.01) ) 
table(FilterGenes_spe)
# 找到满足上面条件的基因：
trait_hubGenes_spe <- colnames(datExpr)[FilterGenes_spe]

write.csv(trait_hubGenes_spe,"trait_hubGenes_spe-04-08-001-93.csv")

pp=automaticNetworkScreeningGS(
  datExpr, GS1, 
  power = 2, networkType = "unsigned", 
  detectCutHeight = 0.995, minModuleSize = min(30, ncol(as.matrix(datExpr))/2), 
  datME = NULL)
write.csv(pp,"pp.csv")
qq=networkScreeningGS(
  datExpr,
  MEs,
  GS1,
  oddPower = 3,
  blockSize = 1000,
  minimumSampleSize =30, 
  addGS = TRUE)
write.csv(qq,"qq.csv")
#Freq_MS_max_cor <- which.max(abs(cor_ADR["Freq",-which(colnames(cor_ADR) == "MEgrey")]))
Freq_MS_max_p <- which.min(p.values["Freq",-which(colnames(p.values) == "MEgrey")])
	

	
###############################

	
	
## step 4 ########################################################
#######################      one step  fa
  
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
	
  table(mergedColors)
  #可以看到每个模块多少个基因
  moduleColors=mergedColors
  # Plot the dendrogram and the module colors underneath
  pdf("step4-genes-modules.pdf")
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  ## assign all of the gene to their corresponding module 
  ## hclust for the genes.



  #明确样本数和基因
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  #首先针对样本做个系统聚类
  datExpr_tree<-hclust(dist(datExpr), method = "average")
  #par(mar = c(0,5,2,0))
  #plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
     #  cex.axis = 1, cex.main = 1,cex.lab=1)
  
  ## 如果这个时候样本是有性状，或者临床表型的，可以加进去看看是否聚类合理
  #针对前面构造的样品矩阵添加对应颜色
 #sample_colors <- numbers2colors(as.numeric(factor(datTraits$ratio)), 
                            #   colors = c("blue","red","green","white","pink"),signed = FALSE)
 sample_colors <- numbers2colors(as.numeric(factor(datTraits$ratio)), 
                                 colors = c("white","red"),signed = FALSE)
  ## 这个给样品添加对应颜色的代码需要自行修改以适应自己的数据分析项目
#sample_colors <- numbers2colors( datTraits ,signed = FALSE)
  ## 如果样品有多种分类情况，而且 datTraits 里面都是分类信息，那么可以直接用上面代码，当然，这样给的颜色不明显，意义不大
  #10个样品的系统聚类树及性状热图
 colnames(datTraits)
 >  colnames(datTraits)
 [1] "D_PVAT"  "ND_PVAT" "PVAT"    "OVF"     "SAF"
 
  par(mar = c(1,4,3,1),cex=0.8)
  #groupLabels =colnames(sample)
  pdf("sample-subtype-cluster.pdf",width=100, height=60, units='mm',res=1200,compression='lzw',pointsize=6)
  plotDendroAndColors(datExpr_tree, sample_colors,
                      groupLabels = "D_PVAT",
                      cex.dendroLabels = 0.5,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.1,
                      main = "Sample dendrogram and trait heatmap")
  dev.off()

  
 ############## 
## step 5  gai gai gai gai   $houmiande 
  #############
## 这一步主要是针对于连续变量，如果是分类变量，需要转换成连续变量方可使用

  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
    design=model.matrix(~0+ datTraits$TAA)
  colnames(design)=levels(datTraits$TAA)
  moduleColors <- labels2colors(net$colors)
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩 (样本vs模块)
  moduleTraitCor = cor(MEs, design , use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  pdf("step5-Module-trait-relationships.pdf",width=60, height=100, units='mm',res=1200,compression='lzw',pointsize=6)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
 # Length of 'xLabels' must equal the number of columns in 'Matrix.'
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()


###############
## step 6
#################
if(T){
  # names (colors) of the modules
  MEs=mergedMEs
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  ## 算出每个模块跟基因的皮尔森相关系数矩
  ## MEs是每个模块在每个样本里面的
  ## datExpr是每个基因在每个样本的表达量
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  modNames
  [1] "blue"        "green"       "black"       "greenyellow" "tan"         "turquoise"   "yellow"     
  [8] "magenta"     "pink"        "purple"      "salmon"      "cyan"        "brown"       "red"        
  [15] "grey" 
  
  ## 只有连续型性状才能只有计算
  ## 这里把是否属 TAA 表型这个变量0,1进行数值化
  #TAA = as.data.frame(design[,2])
  TAA =  datTraits[,1]
  names(TAA) ="TAA"
  
  geneTraitSignificance = as.data.frame(cor(datExpr, TAA , use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.","TAA", sep="");
  names(GSPvalue) = paste("p.GS.", "TAA", sep="");
  
  
  
  module = "blue"
  module = "brown"
  module ="turquoise"
  moduleColors= mergedColors
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  pdf("step6-Module_membership-gene_significance-of-MEgrey.pdf")
  #sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for TAA",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  dev.off()
  
}


## step 7 
if(T){
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  geneTree = net$dendrograms[[1]]; 
  #dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6); 
  dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 20); 
  
  
  

  plotTOM = dissTOM^7; 
  diag(plotTOM) = NA;
  
  
  
  
  
  #TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  nSelect = 4000
  # For reproducibility, we set the random seed
  set.seed(10);
  select = sample(nGenes, size = nSelect);
  selectTOM = dissTOM[select, select];
  # There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
  selectTree = hclust(as.dist(selectTOM), method = "average")
  selectColors = moduleColors[select];
  # Open a graphical window
  sizeGrWindow(9,9)
  # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
  # the color palette; setting the diagonal to NA also improves the clarity of the plot
  plotDiss = selectTOM^7;
  diag(plotDiss) = NA;
  
  png("step7-Network-heatmap4000.png")
  TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
  dev.off()
  
  
  
  #######################################
  # Recalculate module eigengenes
  MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
  
  #
  datKME=signedKME(datExpr, MEs, outputColumnName="kME_MM.")							
  write.csv(datKME, "kME_MM_test.csv")
  
  ## 只有连续型性状才能只有计算
  ## 这里把是否属 TAA 表型这个变量0,1进行数值化
  TAA = as.data.frame(design[,2]);
  names(TAA) = "TAA"
  # Add the weight to existing module eigengenes
  MET = orderMEs(cbind(MEs, TAA))
  ########
  MET = orderMEs(MEs)
  # Plot the relationships among the eigengenes and the trait
  sizeGrWindow(5,7.5);
  
  par(cex = 0.9)
  pdf("step7-Eigengene-dendrogram.pdf")
  plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                        = 90)
  dev.off()
  
  # Plot the dendrogram
  sizeGrWindow(6,6);
  par(cex = 1.0)
  ## 模块的进化树
  pdf("step7-Eigengene-dendrogram-hclust.pdf")
  plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                        plotHeatmaps = FALSE)
  dev.off()
  # Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
  par(cex = 1.0)
  ## 性状与模块热
  
  pdf("step7-Eigengene-adjacency-heatmap.pdf")
  plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                        plotDendrograms = FALSE, xLabelsAngle = 90)
  dev.off()
  
}
  modNames
## step 8 
if(T){
  # Select module
  module ="purple" ;
  module = "turquoise"
  module ="greenyellow"
  
  module = "blue"
  module = "brown"
  module ="grey60"
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  
  gene= modProbes
  write.csv(gene,"grey60.csv")
  write.csv(gene,"turquoise.csv")
  name=modNames
  
  modNames = substring(names(MEs), 3)
  len=length(modNames)

  probes = colnames(datExpr)
  for (j in 1:len){
    module =modNames[j]
    #print( module)
     ## 我们例子里面的probe就是基因
   inModule = (moduleColors==module);
    modProbes = probes[inModule]; 
    gene= modProbes
   fn=paste0('genelist_', module,'.csv')
   write.csv( gene,fn)
  }

  table(net$colors) 
  moduleColors <- labels2colors(net$colors)
  ####################################yes
  a=table( moduleColors)
  write.csv(a,"jishu.csv")
  
  save(TOM,file = 'GSE9106-wgcna-TOM.RData')
  save(probes,moduleColors,net,file = 'GSE9106-wgcna-gene-mode.RData')
## step 9  
if(T){
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(datExpr, power = 20)
  
  
  TOM = TOMsimilarityFromExpr(datExpr, power = 6); 
  # Select module
  module = "brown";
  module = turquoise
  module = "blue"
  module = "brown"
  module ="grey"
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  ## 也是提取指定模块的基因名
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  ## 模块对应的基因关系矩
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  );
}







