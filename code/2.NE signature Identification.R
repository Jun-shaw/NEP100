######################bulk DEGs####################
library(dplyr)
library(ImageGP)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)

matrix <- read.csv('PCaProfilter_matrix.csv')
matrix <- matrix %>%
  group_by(across(1)) %>%  # 根据第一列进行分组
  summarize(across(everything(), max, .names = "{col}"))  # 计算每列的最大值
matrix <- as.data.frame(matrix)
row.names(matrix) <- matrix[,1]
matrix <-matrix[,-1]

group <- read.csv('clinical.csv')
#crpc vs cspc
table(group$PC.type)

crpc_row_names <- group$ID[group$PC.type=='CRPC']
cspc_row_names <- group$ID[group$PC.type=='PRIMARY']

CRPC.matrix <- matrix[,crpc_row_names]
CSPC.matrix <- matrix[,cspc_row_names]

expr <- cbind(CRPC.matrix,CSPC.matrix)
write.csv(expr,'CRPC_vs_primary_matrix.csv')


#差异分析
EZGen:filter_DEGs(expr.data ='CRPC_vs_primary_matrix.csv',
                  data.type ='mRNA',
                  deg.method='DESeq2',
                  Pvalue    =0.05,
                  log2FC    =2,
                  TCGA      =F,
                  tumor.num =484,
                  normal.num=708,
                  color1    ='#CE5C69',
                  color2    ='#5770A6',
                  title     ='PCaProfilter Bulk RNA-seq CRPC vs primary',
                  key.genes =NA)

#对角线出图

CRPC.matrix$CRPC <- rowMeans(CRPC.matrix, na.rm = TRUE) 
CSPC.matrix$Primary <- rowMeans(CSPC.matrix, na.rm = TRUE) 

CRPC.matrix <- tibble::rownames_to_column(CRPC.matrix, var = "Gene_symbol")
CSPC.matrix <- tibble::rownames_to_column(CSPC.matrix, var = "Gene_symbol")

merged_result <- my_result %>%
  left_join(select(CRPC.matrix, Gene_symbol, CRPC), by = "Gene_symbol")
merged_result <- merged_result %>%
  left_join(select(CSPC.matrix, Gene_symbol, Primary), by = "Gene_symbol")
write.csv(merged_result,'mRNA_deseq2.csv')

diffexpr <- read.csv("mRNA_deseq2.CSV",header = T)
diffexpr$CRPC<- log2(diffexpr$CRPC+1)#对平均值进行标准化
diffexpr$Primary<- log2(diffexpr$Primary+1)#对平均值进行标准化
diffexpr$level <- ifelse(diffexpr$padj<0.05, 
                         ifelse(diffexpr$log2FoldChange>=1, "Up", 
                                ifelse(diffexpr$log2FoldChange<=-1, "Down", "NoSig")),"NoSig")#标记差异基因
head(diffexpr)
p <- sp_scatterplot(diffexpr, xvariable = "CRPC", yvariable = "Primary", #定义横纵坐标变量，用的是前面计算的样本平均值
                    color_variable = "level",#颜色以分组定义
                    title ="PCaProfilter Bulk RNA-seq CRPC vs primary", #标题
                    color_variable_order = c("NoSig","Up", "Down"),
                    manual_color_vector = c("grey","#CE5C69","#5770A6")) + #颜色定义
  coord_fixed(1)+ labs(x = "CRPC", y = "Primary")+ theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size = 10))
p
diffexpr$label =""
diffexpr <- diffexpr[order(diffexpr$padj),]
up.genes <- head(diffexpr$Gene_symbol[which(diffexpr$level=="Up")],3)
down.genes <- head(diffexpr$Gene_symbol[which(diffexpr$level=="Down")],3)
top10genes <- c(as.character(up.genes), as.character(down.genes))
diffexpr$label[match(top10genes,diffexpr$Gene_symbol)] <- top10genes
p + geom_text_repel(data = diffexpr, aes(label = label), color = "black", size = 4, fontface = "italic",
                    box.padding = 0.5, point.padding = 0.5, segment.color = 'black', segment.size = 0.3,
                    force = 1, max.iter = 1, min.segment.length = 0.5)





#NEPC vs ARPC
NEPC_row_names <- group$ID[group$Molecular.Type=='NEPC']
ARPC_row_names <- group$ID[group$Molecular.Type=='ARPC']

NEPC.matrix <- matrix[,NEPC_row_names]
ARPC.matrix <- matrix[,ARPC_row_names]

expr <- cbind(NEPC.matrix,ARPC.matrix)
write.csv(expr,'NEPC_vs_ARPC_matrix.csv')

#差异分析
EZGen:filter_DEGs(expr.data ='NEPC_vs_ARPC_matrix.csv',
                  data.type ='mRNA',
                  deg.method='DESeq2',
                  Pvalue    =0.05,
                  log2FC    =1,
                  TCGA      =F,
                  tumor.num =34,
                  normal.num=428,
                  color1    ='#CE5C69',
                  color2    ='#5770A6',
                  title     ='PCaProfilter Bulk RNA-seq NEPC vs ARPC',
                  key.genes =NA)

#对角线出图

NEPC.matrix$NEPC <- rowMeans(NEPC.matrix, na.rm = TRUE) 
ARPC.matrix$ARPC <- rowMeans(ARPC.matrix, na.rm = TRUE) 

NEPC.matrix <- tibble::rownames_to_column(NEPC.matrix, var = "Gene_symbol")
ARPC.matrix <- tibble::rownames_to_column(ARPC.matrix, var = "Gene_symbol")
my_result <- read.csv("mRNA_deseq2.CSV",header = T)

merged_result <- my_result %>%
  left_join(select(NEPC.matrix, Gene_symbol, NEPC), by = "Gene_symbol")
merged_result <- merged_result %>%
  left_join(select(ARPC.matrix, Gene_symbol, ARPC), by = "Gene_symbol")

write.csv(merged_result,'mRNA_deseq2.csv')

diffexpr <- read.csv("mRNA_deseq2.CSV",header = T)
diffexpr$NEPC<- log2(diffexpr$NEPC+1)#对平均值进行标准化
diffexpr$ARPC<- log2(diffexpr$ARPC+1)#对平均值进行标准化
diffexpr$level <- ifelse(diffexpr$padj<0.05, 
                         ifelse(diffexpr$log2FoldChange>=1, "Up", 
                                ifelse(diffexpr$log2FoldChange<=-1, "Down", "NoSig")),"NoSig")#标记差异基因
head(diffexpr)
p <- sp_scatterplot(diffexpr, xvariable = "NEPC", yvariable = "ARPC", #定义横纵坐标变量，用的是前面计算的样本平均值
                    color_variable = "level",#颜色以分组定义
                    title ="PCaProfilter Bulk RNA-seq NEPC vs ARPC", #标题
                    color_variable_order = c("NoSig","Up", "Down"),
                    manual_color_vector = c("grey","#CE5C69","#5770A6")) + #颜色定义
  coord_fixed(1)+ labs(x = "ARPC", y = "ARPC")+ theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size = 10))
p
diffexpr$label =""
diffexpr <- diffexpr[order(diffexpr$padj),]
up.genes <- c('SYP','CHGA','AMIGO2')
down.genes <- c('AR','KLK3','NKX3-1')
top10genes <- c(as.character(up.genes), as.character(down.genes))
diffexpr$label[match(top10genes,diffexpr$Gene_symbol)] <- top10genes
p + geom_text_repel(data = diffexpr, aes(label = label), color = "black", size = 4, fontface = "italic",
                    box.padding = 0.5, point.padding = 0.5, segment.color = 'black', segment.size = 0.3,
                    force = 1, max.iter = 1, min.segment.length = 0.5)


######################WGCNA#######################################
library(WGCNA)
library(tidyverse)
library(tinyarray)
library(sva)
library(limma)

#基础设置
rm(list=ls(all=TRUE))
options(stringsAsFactors = F)# 在读入数据时，遇到字符串后，将其转换成因子，连续型变量要改为FALSE
enableWGCNAThreads()#多线程

# 加载输入文件 
mRNA_TPM_new <-expr
mRNA_TPM_new <- read.csv("PCaProfilter_matrix.csv",row.names = 1)

# 去批次效应
group <- read.csv('clinical.csv')
group <- subset(group, PC.type %in% c('CRPC'))

group$type <-  ifelse(group$Molecular.Type == "NEPC", "NE", "Luminal")

CRPC.matrix <- CRPC.matrix[rowMeans(CRPC.matrix) > 1, ]

CRPC.matrix <- log2(CRPC.matrix+1)
p1 <- draw_pca(exp = CRPC.matrix, group_list = factor(group$Dataset))
p2 <- draw_pca(exp = CRPC.matrix, group_list = factor(group$type))

model <- model.matrix(~factor(group$type))
expr_combat <- ComBat(dat  = CRPC.matrix, batch = group$Dataset,mod = model)

p3 <- draw_pca(exp = expr_combat, group_list = factor(group$Dataset))
p4 <- draw_pca(exp = expr_combat, group_list = factor(group$type))

p1|p2
p3|p4
expr_combat2 <- as.data.frame(expr_combat)

boxplot(expr_combat2[,1:20],main="Non-normalized",las=2)
expr_use <- normalizeBetweenArrays(expr_combat2, method="quantile")
boxplot(expr_use[,1:20],main="Normalized",las=2)

# 开始WGCNA分析
#筛选关键基因
m.mad <- apply(expr_use,1,mad)
expr_use <- expr_use[which(m.mad >max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

#① 对表达矩阵进行转置，行名为样本名，列名为基因
datExpr <- t(expr_use)
#② 判断数据质量，缺失值
gsg <- goodSamplesGenes(datExpr)
gsg$allOK

#③ 判断数据质量，离群样本
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
pdf(file="Figure_1_样本聚类图.pdf",width = 20,height = 4)
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", 
     xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 1.5)
abline(h = 400, col = "#CE5C69") #根据实际情况而定
dev.off()


# 如果有离群值，可以通过以下代码剔除 
clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 30)#自定义
table(clust)
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
datExpr

# 软阈值与平均连接度 
powers = c(seq(1,10,by=1),seq(12,20,by=2))#自定义参数
#调用网络拓扑分析函数
type <- "signed"
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType=type, verbose = 5)
#无尺度的拓扑拟合指数是软阈值功率的一个函数
pdf(file="Figure_2_软阈值2.pdf",width = 10,height = 8)
par(mfrow = c(1,2))
cex1 <-1.0
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="#5770A6");

# 这条线对应的是使用R^2的截止值h
abline(h=0.83,col="#CE5C69",lty=1,lwd=1.5)

# 平均连通性是软阈值功率的一个函数
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="#5770A6")
dev.off()

# 一步法构建共表达网络模型 
power = sft$powerEstimate# 选择一个合适的软阈值，\
nSamples <- ncol(datExpr)
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))
                 )
  )
}
for(i in 1:numCols) {
  datExpr[, i] <- as.numeric(datExpr[, i])
}
net = blockwiseModules(datExpr, #指定用于分析的数据集
                       power = 5,#表示基因之间相关性的绝对值将提高到的幂。更高的幂值将导致网络中基因之间的边的定义更加严格。
                       TOMType = type, #指定要使用的TOM类型。值为“unsigned”，这意味着相关性的绝对值将用于计算TOM。
                       minModuleSize = 200,#指定了一个模块(一组表达高度相关的基因)中必须包含的最小基因数量，这样模块才被认为是重要的。
                       maxBlockSize = 5000,
                       reassignThreshold = 0, #指定将基因重新分配到不同模块所需的模块成员稳定性的最小增量。值为0意味着基因永远不会被重新分配到不同的模块。
                       mergeCutHeight = 0.25,#指定函数生成的树状图(tree-like diagram)被切割以确定模块的高度。较高的值将导致更少、更大的模块。
                       numericLabels = TRUE, #指定输出中的基因标签应该是数字(TRUE)还是字符(FALSE)。
                       pamRespectsDendro = FALSE,#指定围绕中间体(PAM)聚类算法的分区是否应该尊重函数生成的树状图。如果为TRUE, PAM集群将被限制在树状图的分支上。
                       saveTOMs = F,#指定是否将TOM矩阵保存到文件中。
                       saveTOMFileBase = "TPM-TOM", #用于指定保存TOM矩阵时使用的基文件名。
                       verbose = 3)#用于指定函数的详细级别。值越高，输出越详细。

table(net$colors)
mergedColors = net$colors
# 样本的临床数据 
clinical_data <- read.csv("WCGNA_clinical.csv")#读取临床数据
rownames(clinical_data) <- clinical_data$ID#将样本ID作为行名
clinical_data <- clinical_data[,-1]#去除无效信息列
clinical_data_f = clinical_data[rownames(datExpr),]#只提取表达矩阵中样本的临床信息
identical(rownames(datExpr), rownames(clinical_data_f))#判断两个数据行名是否相同

# 识别关键模块 
# 定义基因和样本的数量
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
# 用颜色标签重新计算MEs(计算给定的单个数据集中模块的特征元素（第一主成分）
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes #列名模块颜色
MEs_number <- net$MEs
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitCor#查看模块和临床信息的相关性
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue#计算给定相关性的t test渐近p值。
#作图准备，显示相关性及其p值
sizeGrWindow(16,9)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file="WGCNA.pdf",width = 8,height = 12)
par(mar = c(5, 8, 4, 1))#Bottom margin、Left margin、Top margin、Right margin
#labeledHeatmap用于创建带有标记的行和列的热图
myColors <- colorRampPalette(c('#5770A6', 'white', '#CE5C69'))(100)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),#临床信息文件列名
               yLabels = names(MEs),
               ySymbols = names(MEs),#行的标签
               colorLabels = FALSE,#指示是否显示颜色标签的逻辑值
               colors = myColors,#指定热图中使用的颜色
               border_color = "write",
               textMatrix = textMatrix,#要添加到热图中的文本标签矩阵
               setStdMargins = FALSE,#指示是否使用标准边距的逻辑值
               cex.text = 1,#指定文本标签的大小
               zlim = c(-1,1),#指定要在热图中绘制的值的范围
               main = paste("Module-trait relationships"),#主标题
               xLabelsAngle = 45,#指定绘制x轴标签的角度
               yLabelsPosition = 'left',#指定y轴标签位置的字符串(例如，'left'或'right')
               x.adj.lab.y = 1#指定x轴标签的垂直位置
)
dev.off()
up.module <- c(3,7,8)
down.module <- c(2,5)
up_WGCNA <- colnames(datExpr)[moduleColors%in%up.module]
down_WGCNA<- colnames(datExpr)[moduleColors%in%down.module]

write.csv(up_WGCNA,'WCGNAupgenes.csv')
write.csv(down_WGCNA,'WCGNAupgenes.csv')

######################findmarker##################
seurat.harmony <- readRDS(seurat.harmony)

NE.markers <- FindMarkers(seurat.harmony, ident.1 = "NE", ident.2 = "Luminal",logfc.threshold = 0.5)
write.csv(NE.markers,'FindMarkers_DEG(NE VS L).csv')

NE.markers2 <- FindMarkers(seurat.harmony, ident.1 = "NE",logfc.threshold = 0.5)
write.csv(NE.markers2,'FindMarkers_DEG(NE VS ALL).csv')


######################伪bulk###################
bs = split(colnames(seurat.harmony),seurat.harmony@meta.data$sampleID)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(seurat.harmony) %in% bs[[x]]
    rowSums(as.matrix(seurat.harmony@assays$RNA@layers$counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
rownames(ct) <- rownames(seurat.harmony)
# 获取分组信息
phe = seurat.harmony@meta.data[,c('sampleID','NE')]
phe = unique(seurat.harmony@meta.data[,c('sampleID','NE')]) 
phe[1:5,1:2]
group_list = phe[match(names(bs),phe$sampleID),'NE']
group_list <- factor(group_list,levels = c("NE","nonNE"))
table(group_list) 

# 赋值并对每一行的
exprSet = ct
dim(exprSet)
exprSet=exprSet[apply(exprSet,1, function(x) sum(x>1) > 1),]

# DESeq2分析
colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
# exprSet <- apply(exprSet, 2, as.integer) 这一步不需要哦，可以想一想
# rownames(exprSet) <- rownames(scRNA) 如果这步上面没有做，可以再加上去
dds <- DESeqDataSetFromMatrix(countData = exprSet,
                              colData = colData,
                              design = ~ group_list)
dds <- DESeq(dds) 

res <- results(dds, 
               contrast=c("group_list",
                          levels(group_list)[1],
                          levels(group_list)[2]))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG_deseq2 = na.omit(DEG)

write.csv(DEG_deseq2,'Pseudobulk_DEGs.csv')
######################vene#################################

CRPC_Primary <- read.csv('CRPC vs primary/mRNA_DEGs.csv')
up_CRPC_Primary <- CRPC_Primary$Gene_symbol[CRPC_Primary$regulate=='up-regulated']
down_CRPC_Primary <- CRPC_Primary$Gene_symbol[CRPC_Primary$regulate=='down-regulated']

NEPC_ARPC <- read.csv('NEPC vs ARPC/mRNA_DEGs.csv')

up_NEPC_ARPC <- NEPC_ARPC$Gene_symbol[NEPC_ARPC$regulate=='up-regulated']
down_NEPC_ARPC <- NEPC_ARPC$Gene_symbol[NEPC_ARPC$regulate=='down-regulated']

up_WGCNA <- read.csv('WCGNAupgenes.csv')[,1]
down_WGCNA <- read.csv('WCGNAdowngenes.csv')[,1]

FindMarkers <- read.csv('G:/importance/undergraduated/2_CRPC/第二版/单细胞测序建模/FindMarkers_DEG(NE VS ALL).csv')
up_FindMarkers<- FindMarkers$X[FindMarkers$avg_log2FC >= 2&FindMarkers$p_val_adj <=0.01&FindMarkers$PCTD>=0.2]
FindMarkers2 <- read.csv('G:/importance/undergraduated/2_CRPC/第二版/单细胞测序建模/FindMarkers_DEG(NE VS L).csv')
down_FindMarkers<- FindMarkers2$X[FindMarkers2$avg_log2FC <= -2&FindMarkers2$p_val_adj <=0.01&FindMarkers2$PCTD>=0.25]

Pseudobulk <- read.csv('G:/importance/undergraduated/2_CRPC/第二版/单细胞测序建模/Pseudobulk_DEGs.csv')
up_Pseudobulk<- Pseudobulk$X[Pseudobulk$log2FoldChange >= 2&Pseudobulk$padj <=0.01]
down_Pseudobulk<- Pseudobulk$X[Pseudobulk$log2FoldChange <= -2&Pseudobulk$padj <=0.01]

library (VennDiagram)

up_set1=t(up_WGCNA)
up_set2=t(up_CRPC_Primary)
up_set3=t(up_FindMarkers)
up_set4=t(up_Pseudobulk)
up_set5=t(up_NEPC_ARPC)

venn.diagram(x=list(up_set1,up_set2,up_set3,up_set4,up_set5),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('white','white','white','white','white'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex =1.2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#CE5C69','#F4C889','#BDD5A3','#A281B1','#8BA1C6'), # 填充色 配色https://www.58pic.com/
             category.names = c("WGCNA\nNE related", 'Bulk data\nCRPC vs CSPC','FindMarkers\nNE related','Pseudobulks\nNEPC related','Bulk data\n NEPC vs ARPC') , #标签名
             cat.dist = 0.08, # 标签距离圆圈的远近
             cat.pos = c(-0, -0,-180, -180,-0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex =1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "text",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='G:/importance/undergraduated/2_CRPC/第二版/DEG、WGNA/upgene_VENE.TIF',# 文件保存
             imagetype="tiff",  # 类型（tiff png svg）
             resolution = 600,  # 分辨率
             compression = "lzw"# 压缩算法
)
intersection_up_sample <- intersect(up_set1, up_set2)
intersection_up_sample <- intersect(intersection_up_sample, up_set3)
intersection_up_sample <- intersect(intersection_up_sample, up_set4)
NEPupsig <- intersect(intersection_up_sample, up_set5)
intersection_up_sample
write.csv(NEPupsig,'NEPupsig.csv')

down_set1=t(down_WGCNA)
down_set2=t(down_CRPC_Primary)
down_set3=t(down_FindMarkers)
down_set4=t(down_Pseudobulk)
down_set5=t(down_NEPC_ARPC)

venn.diagram(x=list(down_set1,down_set2,down_set3,down_set4,down_set5),
             scaled = F, # 根据比例显示大小
             alpha= 0.5, #透明度
             lwd=1,lty=1,col=c('white','white','white','white','white'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             cex = 1.2, # 数字大小
             fontface = "bold",  # 字体粗细；加粗bold
             fill=c('#B30C2A','#E0A980','#519981','#735c88','#5770A6'), # 填充色 配色https://www.58pic.com/
             category.names = c("WGCNA\nLunimal related", 'Bulk data\nCSPC vs CRPC ','FindMarkers\nLunimal related','Pseudobulks\nCSPC related','Bulk data\nARPC vs NEPC') , #标签名
             cat.dist = 0.08, # 标签距离圆圈的远近
             cat.pos = c(-0, -0,-180, -180,-0), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             cat.cex =1, #标签字体大小
             cat.fontface = "bold",  # 标签字体加粗
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             cat.default.pos = "text",  # 标签位置, outer内;text 外
             output=TRUE,
             filename='G:/importance/undergraduated/2_CRPC/第二版/DEG、WGNA/downgene_VENE.TIF',# 文件保存
             imagetype="tiff",  # 类型（tiff png svg）
             resolution = 600,  # 分辨率
             compression = "lzw"# 压缩算法
)


intersection_down_sample <- intersect(down_set1, down_set2)
intersection_down_sample <- intersect(intersection_down_sample, down_set3)
intersection_down_sample <- intersect(intersection_down_sample, down_set4)
intersection_down_sample <- intersect(intersection_down_sample, down_set5)
NEPdownsig
write.csv(NEPdownsig,'NEPdownsig.csv')


######################计算基因集评分###########
genesets <- read.csv('G:/importance/undergraduated/数据集整理/genesets.csv')
saveRDS(genelist,'genelist.rds')
Cal_score <- function(seurat.harmony=seurat.harmony,genesets=genesets,ctrl = 100,seed=520){
  
  genelist_all <- list()
  for (i in seq_along(genesets)) {
    library(Seurat)
    genelist.name <- colnames(genesets)[i]
    print(genelist.name)
    genelist<- genesets[,i]
    genelist <- genelist[nzchar(genelist)]
    
    genelist_all[[length(genelist_all) + 1]] <- genelist
    names(genelist_all)[i] <- genelist.name
    
    seurat.harmony <- AddModuleScore(seurat.harmony,
                                     features = list(genelist),
                                     ctrl = ctrl,search = T,seed = seed,
                                     name = genelist.name)
    
    colname <- grep(genelist.name, colnames(seurat.harmony@meta.data), value = TRUE)
    if (length(colname) > 1) {
      seurat.harmony@meta.data[[genelist.name]] <- rowMeans(seurat.harmony@meta.data[, colname], na.rm = TRUE)
      seurat.harmony@meta.data[, colname] <- NULL
    } else if (length(colname) == 1) {
      seurat.harmony@meta.data[[genelist.name]] <- seurat.harmony@meta.data[[colname]]
      seurat.harmony@meta.data[[colname]] <- NULL
    }
    
    #parts <- strsplit(genelist.name, "\\.")[[1]]
    #genelist.name_up_down <- paste(parts[1:2], collapse = ".")
    #col_up_down <- grep(genelist.name_up_down, colnames(seurat.harmony@meta.data), value = TRUE)
    #if (length(col_up_down) == 2) {
    #  seurat.harmony@meta.data[[genelist.name_up_down]] <- seurat.harmony@meta.data[,col_up_down[[1]]]-seurat.harmony@meta.data[,col_up_down[[2]]]
    #  seurat.harmony@meta.data[, col_up_down] <- NULL
    #}
    
  }
  
  return(seurat.harmony@meta.data)
}

score.res <- Cal_score(seurat.harmony=seurat.harmony,genesets=genesets)
score.res <- seurat.harmony@meta.data

#####################计算celltype的平均值
celltype_values <- unique(score.res$celltype12)
result_matrix <- matrix(NA, nrow = length(celltype_values), ncol =29)
rownames(result_matrix) <- celltype_values
for (i in 1:29) {
  for (celltype in celltype_values) {
    subset_data <- score.res[score.res$celltype12 == celltype, i+30]
    result_matrix[celltype, i] <- mean(subset_data, na.rm = TRUE)
  }
}
colnames(result_matrix) <- colnames(score.res)[31:59]

View(result_matrix)
result_matrix <- t(result_matrix)
write.csv(result_matrix,'celltype_score.csv')
###################计算sampletype的平均值
sampleID_values <- unique(score.res$sampleID)
result_matrix2 <- matrix(NA, nrow = length(sampleID_values), ncol =29)
rownames(result_matrix2) <- sampleID_values
for (i in 1:29) {
  for (sampleID in sampleID_values) {
    subset_data <- score.res[score.res$sampleID == sampleID, i+30]
    result_matrix2[sampleID, i] <- mean(subset_data, na.rm = TRUE)
  }
}
colnames(result_matrix2) <- colnames(score.res)[31:59]

View(result_matrix2)
result_matrix2 <- t(result_matrix2)
write.csv(result_matrix2,'sampleID_score.csv')
######################小提琴图##################
library(ggpubr) 
for (i in seq_along(genelist_all)) {
  print(i)
  geneSet <- names(genelist_all)[i]
ggviolin(seurat.harmony@meta.data, x="celltype", y=geneSet, width = 1, 
         color = "black",#轮廓颜色
         fill="celltype",#填充
         palette = col,
         add = 'mean_sd',
         ylab=geneSet,font.y = 24,
         xlab = F, #不显示x轴的标签
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size=1, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "none"
)+ theme(
  axis.text.x = element_text(size = 20, angle = 45, hjust = 1),  # x轴标签设置
  axis.text.y = element_text(size = 20))  # y轴标签设置
ggsave(filename = paste0('AUCell/', geneSet, '.pdf'), width = 6, height = 6, units = 'in')
}
######################散点图####################
library(tidyverse)
library(ggplot2)
library(ggpubr)
data <- read.csv('sampleID_score.csv',row.names = 1)
colnames(data)
data$sampletype <- factor(data$sampletype,levels = c('Primary','mCSPC','CRPC','mCRPC','NEPC'))

datane <- data %>%
  filter(NE == 'NE')
# 根据 NE 的不同值指定颜色
for (i in 1:27) {
  datasetID <- colnames(data)[i]
  print(datasetID)
  #############散点图1
  ggplot() +
    geom_point(data = data,
               aes(x = ISUP, y = data[[datasetID]], fill = NE),
               size = 5,
               shape = 21,
               color = "black") +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.text = element_text(face = "italic",size = 24),
      axis.title.x = element_text(size = 24), axis.text = element_text(size = 20,color = "black"),
      axis.title.y = element_text(size = 24),axis.line = element_line(color = "black", size = 0),
      panel.border = element_rect(color = "black", fill = NA, size = 2)  # 添加边框
    ) +
    guides(fill = guide_legend(ncol = 2, title = NULL,size = 24)) +
    labs(x = "ISUP",
         y = "Mudule Score") +
    scale_fill_manual(values = c("NE" = "#ce5c69", "nonNE" = "#5770a6")) +
    stat_cor(data = data, aes(x = ISUP, y = data[[datasetID]]), method = "pearson", 
             label.y = max(data[[datasetID]], na.rm = TRUE) * 0.9, size = 6) +
    # 添加带置信区间的回归线
    geom_smooth(data = data, aes(x = ISUP, y = data[[datasetID]]), 
                method = "lm", 
                se = TRUE,  # 显示置信区间
                color = "#a281b1", 
                linetype = "dashed",  # 设置为虚线
                size = 1.5)  # 设置线条粗细
  ggsave(filename = paste0('ISUP/', datasetID, '.pdf'), width = 6, height = 6, units = 'in')
}
######################箱式图#####################
for (i in 1:27) {
  datasetID <- colnames(data)[i]
  print(datasetID)
  
  comparisons_list <- list(
    c("NEPC", "CRPC"),
    c("NEPC", "mCRPC"),
    c("NEPC", "mCSPC") ,
    c("NEPC", "Primary") # 添加更多的比较组
  )
  
  # 使用 ggplot 绘制图形
  ggplot(seurat.harmony@meta.data, aes(x = sampletype, y = datasetID)) +
    geom_boxplot(outlier.size = 1.5, size = 1.2, aes(color = NE), outlier.shape = NA) +  # 设置箱式图的边框颜色
    geom_jitter(aes(fill = NE), position = position_jitter(width = 0.3), size = 5, shape = 21, color = "black") +  # 设置散点的边框颜色为黑色
    theme_minimal(base_size = 16, base_family = "sans") + 
    theme(
      text = element_text(size = 16),  # 调整字体大小
      axis.title = element_text(size = 24),  # 调整轴标题的大小
      axis.text = element_text(size = 20,color = "black"), 
      panel.grid.major = element_blank(),  # 去掉主网格线
      panel.grid.minor = element_blank(),  # 去掉次网格线
      axis.line = element_line(size = 0.5)  # 添加坐标轴线
    ) +
    labs(x = "", y = "Module Score") +
    scale_color_manual(values = c("#ce5c69", "#5770a6"), name = NULL) +  # 设置箱式图的边框颜色
    scale_fill_manual(values = c("#ce5c69", "#5770a6"), name = NULL) +  # 设置散点的填充色
    geom_signif(comparisons = comparisons_list, test = "wilcox.test", map_signif_level = TRUE) +  # 添加多个组的比较
    facet_wrap(~ datasetID, scales = "free_y")  # 按 datasetID 分面，y轴自由缩放
  
  ggsave(filename = paste0('sampletype/', datasetID, '.pdf'), width = 6 ,height = 6, units = 'in')
}

######################热图#########
bs = split(colnames(sce),sce@meta.data$celltype12)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(sce) %in% bs[[x]]
    rowSums(as.matrix(sce@assays$RNA@layers$counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
rownames(ct) <- rownames(sce)

z_trans <- function(data=ct){
  # 计算 Z-scores
  z_scores <- t(apply(data, 1, function(x) {
    (x - mean(x)) / sd(x)
  }))
  
  # 将 Z-scores 转换为数据框
  z_scores_df <- as.data.frame(z_scores)
  rownames(z_scores_df) <- rownames(data)
  colnames(z_scores_df) <- colnames(data)
  
  # 归一化到 -2 到 2 之间
  # 计算 Z-score 的最小值和最大值
  z_min <- min(z_scores_df)
  z_max <- max(z_scores_df)
  
  # 线性变换
  normalized_scores <- 4 * ((z_scores_df - z_min) / (z_max - z_min)) - 2
  
  # 将结果转换为数据框
  normalized_scores_df <- as.data.frame(normalized_scores)
  rownames(normalized_scores_df) <- rownames(data)
  colnames(normalized_scores_df) <- colnames(data)
  
  # 结果矩阵
  result_matrix <- normalized_scores_df
  return(result_matrix)
}
result_matrix <- z_trans(ct)

row_names <- rownames(result_matrix)

# 创建一个分组向量，根据行名分为 NEP100.up 和 NEP100.down
row_split <- ifelse(row_names%in%NEP100.up, "NEP100.up", "NEP100.down")

row_split <- factor(row_split,levels = c("NEP100.up", "NEP100.down"))


order_index <- order(result_matrix$NE, decreasing = TRUE)  # 从大到小排序

# 根据排序后的索引重新排列 result_matrix



sorted_result_matrix <- result_matrix[order_index, ]



cluster_colors<-c('#b30c2a','#ce5c69','#F99999','#e0a980','#f4c889','#bdd5a3','#86a979','#519981','#8ba1c6','#5770a6','#a281b1','#735c88')
# 使用 result_matrix 的列名作为注释

ha_top <- HeatmapAnnotation(
  Type = anno_simple(c('Epithelial', 'Luminal', 'NE', 'Fib', 'MyoFib', 
                       'Ecs', 'Mast', 'Mye', 'T', 'NK', 'B', 'Plasma'),
                     col = setNames(cluster_colors, 
                                    c('Epithelial', 'Luminal', 'NE', 'Fib', 'MyoFib', 
                                      'Ecs', 'Mast', 'Mye', 'T', 'NK', 'B', 'Plasma')),  # 映射颜色
                     border = TRUE
  )
)

category_vector <- factor(c(rep('NEP100.up', 90), rep('NEP100.down', 10)),
                          levels = c('NEP100.up', 'NEP100.down'))
ha_left <- rowAnnotation(
  Type = anno_simple(
    category_vector,
    col = setNames(c('#5770a6', '#ce5c69'), 
                   levels(category_vector)),
    border = TRUE
  )
)

plot(ha_left)

celltype <- factor(c('Epithelial', 'Luminal', 'NE', 'Fib', 'MyoFib', 
                     'Ecs', 'Mast', 'Mye', 'T', 'NK', 'B', 'Plasma'),
                   levels = c('Epithelial', 'Luminal', 'NE', 'Fib', 'MyoFib', 
                              'Ecs', 'Mast', 'Mye', 'T', 'NK', 'B', 'Plasma'))

col_fun = colorRamp2(c(-2, 0, 2), c("#364888", "white", "#b30c2a"))
col_fun = colorRamp2(c(-2, 0, 2), c("#5770a6", "white", "#b30c2a"))
p1 <- Heatmap(sorted_result_matrix,  # 不包含 NE 列
              name = "Module Score",
              col = col_fun,
              na_col = "grey",
              border = TRUE,
              border_gp = gpar(lty = 1, lwd = 3, col = "black"),
              # rect_gp = gpar(col = "white"),
              cluster_rows = F,
              cluster_columns = FALSE,
              row_names_side = "left",
              row_title = "",
              row_title_side = "left",
              row_title_gp = gpar(fontsize = 16, fontface = "bold"),
              show_row_names = F,
              column_title = "",
              column_title_side = "top",
              column_title_gp = gpar(fontsize = 14, fontface = "bold"),
              top_annotation = ha_top,
              left_annotation=ha_left,
              column_split = celltype,column_gap = unit(2, "mm"),
              row_split = row_split,row_gap = unit(6, "mm"),
              heatmap_legend_param = list(
                title = "Z-score",  # 图例标题
                title_gp = gpar(fontsize = 18, fontface = "bold"),  # 标题样式
                labels_gp = gpar(fontsize = 14),  # 标签样式
                grid_width = unit(5, "mm"),  # 调整图例的宽度
                grid_height = unit(10, "mm")  # 调整图例的高度
              ))


bs2 = split(colnames(sce),sce@meta.data$sampletype)
ct2 = do.call(
  cbind,lapply(names(bs2), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(sce) %in% bs2[[x]]
    rowSums(as.matrix(sce@assays$RNA@layers$counts[, kp]  ))
  })
)

colnames(ct2) <- names(bs2)
rownames(ct2) <- rownames(sce)

result_matrix2 <- z_trans(ct2)
desired_order <- c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC')
sorted_result_matrix2 <- result_matrix2[rownames(sorted_result_matrix), desired_order]
sampletype_colors <- c('#ce5c69','#f4c889','#bdd5a3','#5770a6','#a281b1')
sampletype <- factor(c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC'),
                     levels = c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC'))
ha_top2 <- HeatmapAnnotation(
  Type = anno_simple(c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC'),
                     col = setNames(sampletype_colors, 
                                    c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC')),  # 映射颜色
                     border = TRUE
  )
)
p2 <- Heatmap(sorted_result_matrix2,  # 不包含 NE 列
              name = "Module Score",
              col = col_fun,
              na_col = "grey",
              border = TRUE,
              border_gp = gpar(lty = 1, lwd = 3, col = "black"),
              # rect_gp = gpar(col = "white"),
              cluster_rows = F,
              cluster_columns = FALSE,
              row_names_side = "left",
              row_title = "",
              row_title_side = "left",
              row_title_gp = gpar(fontsize = 16, fontface = "bold"),
              show_row_names = F,
              column_title = "",
              column_title_side = "top",
              column_title_gp = gpar(fontsize = 14, fontface = "bold"),
              top_annotation = ha_top2,
              column_split = sampletype,
              column_gap = unit(2, "mm"),
              row_split = row_split,
              row_gap = unit(6, "mm"),
              heatmap_legend_param = list(
                title = "Z-score",  # 图例标题
                title_gp = gpar(fontsize = 18, fontface = "bold"),  # 标题样式
                labels_gp = gpar(fontsize = 14),  # 标签样式
                grid_width = unit(5, "mm"),  # 调整图例的宽度
                grid_height = unit(10, "mm")  # 调整图例的高度
              ))
p1
p2