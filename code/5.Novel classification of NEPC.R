#########################################
##################Figure 6A-B##############
#########################################
NE_sce = seurat.harmony[, Idents(seurat.harmony) %in% 'NE']

#标准化等
NE_sce <- scRemoveBatch(seurat.data = NE_sce, n.pcs = 30,methods = "harmony",batchID = "sampleID" )

#选参
ElbowPlot(NE_sce, ndims=30, reduction="harmony") 
pct <- NE_sce[["harmony"]]@stdev / sum(NE_sce[["harmony"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
pc.use <- which(cumu >= 90)[1]
#降维聚类
NE_sce <- FindNeighbors(NE_sce,reduction = 'harmony',dims = 1:pc.use)#选择多少个PCA维度进行分析

NE_sce <- FindClusters(NE_sce,resolution = seq(from = 0.1,to = 0.5, by = 0.1))#测试选择分辨率大小

clustree(NE_sce) 

NE_sce <- FindClusters(NE_sce,resolution = 0.1)

#可视化

NE_sce <- RunUMAP(NE_sce,dims = 1:14,reduction = 'harmony')

#默认为active.ident
DimPlot(NE_sce, reduction = "umap", label = TRUE,group.by = "seurat_clusters", pt.size=1.4)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank())Idents(NE_sce)<- NE_sce$seurat_clusters

######计算NEP100表达

av.exp<- AggregateExpression(NE_sce)$RNA
features=names(tail(sort(apply(av.exp, 1, sd)),200))

av.exp<- av.exp[which(row.names(av.exp)%in% NEP100),]
av.exp <- as.data.frame(av.exp)
av.exp <- cor(av.exp, method= "spearman")
pheatmap::pheatmap(av.exp)

subtype1<- c(0)
subtype2<- c(1,8)
subtype3<- c(2,5,6,7,3)
subtype4<- c(4)

NE_sce@meta.data <- NE_sce@meta.data %>% 
  mutate(subtype = case_when(
    seurat_clusters %in% subtype1 ~ 'subtype1',
    seurat_clusters %in% subtype2 ~ 'subtype2',
    seurat_clusters %in% subtype3 ~ 'subtype3',
    seurat_clusters %in% subtype4 ~ 'subtype4'))

saveRDS(NE_sce,'NE_sce.rds')
#########################################
##################Figure 6C##############
#########################################

library(ggplot2)
library(dplyr)
library(msigdbr)
library(Seurat)
library(GSVA)
library(pheatmap)
library(patchwork)
expr <- GetAssayData(NE_sce) #表达矩阵
meta <- NE_sce@meta.data[,c("seurat_clusters","subtype")] #类别
m_df = msigdbr(species = "Homo sapiens", category = "H")#, subcategory = "CP:KEGG") #选取物种人类
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

pathway <- read.csv('pathway.csv')
genelist <- path_f(pathway)

expr=as.matrix(expr) 

kegg <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10) #gsva
hallmark <- gsva(expr, msigdbr_list, kcdf="Gaussian",method = "gsva",parallel.sz=10)
pathway_DIY<- gsva(expr, genelist, kcdf="Gaussian",method = "gsva",parallel.sz=10)

saveRDS(kegg,'kegg.rds')
saveRDS(hallmark,'hallmark.rds')
saveRDS(NE_sce,'NE_sce.rds')

##每个细胞类别与功能相关热图
meta <- meta %>%arrange(meta$subtype)
data <- hallmark[,rownames(meta)]
group <- factor(meta[,"subtype"],ordered = T)
data1 <-NULL
subtype <- unique(group)
for(i in seq_along(subtype)){
  ind <-which(group==subtype[i])
  dat <- apply(data[,ind], 1, mean)
  data1 <-cbind(data1,dat)
  }

colnames(data1)<-subtype

Zscore <- function(expr=av.exp,feature=feature){
  z_scores <-  t(scale(t(expr)))
  z_scores_normalized <- pmax(pmin(z_scores, 2), -2)
  z_scores_normalized <- as.matrix(z_scores_normalized)
  return(z_scores_normalized)
}

data1 <- Zscore(data1)
new_min <- -2
new_max <- 2

# 计算 z_scores 的最小值和最大值
old_min <- min(data1)
old_max <- max(data1)

# 进行线性变换
z_scores_scaled <- (data1 - old_min) / (old_max - old_min) * (new_max - new_min) + new_min



library(pheatmap)
p <- pheatmap(z_scores_scaled,
              cluster_rows = F,
              cluster_cols = F,
              show_rownames = T,
              show_colnames = T,
              color =colorRampPalette(c("#5770A6", "white","#CE5C69"))(100),
              cellwidth = 15, cellheight = 15,
              fontsize = 10)
pdf(("gsva_celltype.pdf"),width = 8,height = 40)
print(p)
dev.off()
#########################################
##################Figure 6E##############
#########################################
genes <- c(#Neur
  'CHGA','CHGB','STX19','SNAP25','CA9',
  #Prol
  'STMN1','EZH2','MKI67','PCNA','TOP2A',
  #EMT,
  'S100A6','S100A11','TWIST1','VIM','CDH2'
)

DotPlot(object = NE_sce,features = genes,group.by = 'subtype') + coord_flip()


#########################################
##################Figure 6G##############
#########################################
GO <- enrichGO(gene = Genes$ENTREZID, #输入基因的"ENTREZID"
               OrgDb = org.Hs.eg.db,#注释信息
               keyType = "ENTREZID",
               ont = "all",     #可选条目BP/CC/MF
               pAdjustMethod = "BH", #p值的校正方式
               pvalueCutoff = 1,   #pvalue的阈值
               qvalueCutoff = 1, #qvalue的阈值
               minGSSize = 5,
               maxGSSize = 5000,
               readable = TRUE)   #是否将entrez id转换为symbol
GO_result <- as.data.frame(GO)#转化结果
barplot(GO)
write.csv(GO,'subtype1_GO.csv')

#########################################
##################Figure 6H##############
#########################################
library(CaCTS)

Idents(NE_sce)<- NE_sce$seurat_clusters
av.exp<- AggregateExpression(NE_sce)$RNA
av.exp <- as.data.frame(av.exp)

av.exp$gene <- row.names(av.exp)

group<- read.csv("group.csv")
colnames(group)[1] <- "sample.id"
colnames(group)[2] <- "group.name"

TF.list = read.delim("merged.list.1671.TFs.txt", sep = "\t")

f.TCGA.RNA = av.exp[which(av.exp$gene %in% as.character(TF.list$NameTF)),]
dim(f.TCGA.RNA)

aux = which(!as.character(TF.list$NameTF) %in% f.TCGA.RNA$gene)
write.table(TF.list[aux,],file = "non-expressed-TFs.txt", quote = F, row.names = F)

f.TCGA.RNA <- f.TCGA.RNA[,3:ncol(f.TCGA.RNA)]
f.TCGA.RNA <- f.TCGA.RNA[,-10]

delta1 = max(f.TCGA.RNA, na.rm = T) - min(f.TCGA.RNA, na.rm = T)
delta2 = max(f.TCGA.RNA, na.rm = T) - 0

f.TCGA.RNA.rs = (f.TCGA.RNA - min(f.TCGA.RNA, na.rm = T)) * delta1 / delta2
f.TCGA.RNA.rs[1:5, 1:4]


matrix.rep <- prepare_representaive_samples(expr.matrix = f.TCGA.RNA.rs, sample.descr = group, save.file = F)
matrix.rep[1:4, 1:4]

res.CaCTS1 <- run_CaCTS_score(matrix.rep,'CNE')
res.CaCTS2 <- run_CaCTS_score(matrix.rep,'Prol1')
res.CaCTS3 <- run_CaCTS_score(matrix.rep,'Prol2')
res.CaCTS4 <- run_CaCTS_score(matrix.rep,'EMT')

matching_genes <- rownames(NE_sce) %in% res.CaCTS4$Name
NE_genes <- subset(NE_sce, features = rownames(NE_sce)[matching_genes])

filtered1 <- filter_by_expression_rank(rep.matrix = matrix.rep, tf.scores = res.CaCTS, query.name = cancer, pn=0.05, pnE=0.05)
filtered2 <- filter_by_expression_rank(rep.matrix = matrix.rep, tf.scores = res.CaCTS, query.name = cancer, pn=0.05, pnE=0.05)
filtered3 <- filter_by_expression_rank(rep.matrix = matrix.rep, tf.scores = res.CaCTS, query.name = cancer, pn=0.05, pnE=0.05)
filtered4 <- filter_by_expression_rank(rep.matrix = matrix.rep, tf.scores = res.CaCTS, query.name = cancer, pn=0.05, pnE=0.05)

filtered1$subtype <- 'CNE'
filtered2$subtype <- 'Prol1'
filtered3$subtype <- 'Prol2'
filtered4$subtype <- 'EMT'


filtered <- rbind(filtered1,filtered2,filtered3,filtered4)
write.csv(filtered,'filtered.csv')


TF <- read.csv('filtered.csv')


library(scplotter)
FeatureStatPlot(NE_sce, features = Name, ident = "subtype_total2", plot_type = "dot")



ggplot(TF, aes(x = subtype, y = Name, size = LogValue, fill = Expr.mean)) +
  geom_point(shape = 21, color = "black") +  # 使用黑色边框
  scale_size(range = c(1, 20), name = "Log Value") +  # 调整气泡大小范围
  scale_fill_gradient(low = "#5770A6", high = "#CE5C69", name = "Expr Mean") +  # 设置颜色渐变
  theme_minimal() +
  theme(
  ) + 
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +  # 轴标签
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend("Percent Expression")) +  # legend
  scale_color_gradientn(colours = c('#5770a6',"white", "#ce5c69"),limits = c(-1.5, 1.5))+
  guides(size = guide_legend(title = "Percent Expression"))# 颜色   # 去掉# 添加标题

#########################################
##################Figure 6I##############
#########################################

Idents(NE_sce)<- NE_sce$seurat_clusters
av.exp<- AggregateExpression(NE_sce)$RNA
av.exp <- as.data.frame(av.exp)
feature <- c('ONECUT2','HOXB5','NEUROD1','HES6','POU2F3','E2F1','YAP1','ID3',
             'CHGA', 'SYP' ,'ACTL6B', 'SNAP25' ,'INSM1', 'CHRNB2' ,'SRRM4' ,'ASCL1',
             'CELF3' ,'PCSK1' ,'SOX2' ,'POU3F2' ,'LMO3' ,
             'AR' ,'NKX3-1' ,'KLK3', 'SPDEF', 'TMPRSS2',
             'VIM','CDH2','MCAM','FN1')

z_scores_normalized <- Zscore(av.exp,feature)

col_fun = colorRamp2(c(-2,0,2), c('#5770a6','white',"#ce5c69"))

new_order <- feature
new_order2 <- c('g0','g6','g5','g3','g2','g7','g4','g1','g8')
new_order2 <- c('CNE','Prol1','Prol2','EMT')
z_scores_normalized <- z_scores_normalized[match(new_order, row.names(z_scores_normalized)), match(new_order2, colnames(z_scores_normalized))]


row_split <- factor(c(rep('TF',8),rep("NE1", 8), rep("NE2", 5), rep("AR", 5), rep("SQU", 4)),
                    levels = c('TF',"NE1", "NE2", "AR", "SQU"))
col_split <- factor(c("subtype1", rep("subtype2", 5), "subtype3", rep("subtype4", 2)),
                    levels = c("subtype1", "subtype2", "subtype3", "subtype4"))

write.csv(z_scores_normalized,'1.csv')

z_scores_normalized <- read.csv('1.csv',row.names = 1)
z_scores_normalized <- as.matrix(z_scores_normalized)
Heatmap(z_scores_normalized,  # 不包含 NE 列
        name = "Z-score",
        col = col_fun,
        na_col = "grey",
        border = TRUE,
        border_gp = gpar(lty = 1, lwd = 2, col = "black"),
        rect_gp = gpar(col = "white"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "right",
        row_title = "",
        row_split = row_split,
        column_split = col_split,
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_names = TRUE,
        column_title = "",  # 添加列标题
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"))
#row_gap = unit(5, "mm"),  # 调整行分割的宽度
#column_gap = unit(5, "mm"),  # 调整列分割的宽度
# right_annotation = ha_right,
#top_annotation = ha_top_combined)
#########################################
##################Figure J-K##############
#########################################
#see the code’1.Single-cell atlas.R‘
#see the code’1.Single-cell atlas.R‘  

#########################################
##################Figure 7A##############
#########################################
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(RColorBrewer)
library(future)
library(clustree)
library(cowplot)
library(stringr)
library(monocle3)
plan("multisession", workers = 4)
plan()
options(future.globals.maxSize = 10 * 1024^3)

out.prefix = "./Outplot/"

#### 1.input data
NE.data = read_rds("./NE_sce.rds")

#monocle3
### 1.1 Make the CDS object
data <- GetAssayData(NE.data, assay = 'RNA', slot = 'counts')
cell_metadata <- NE.data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

### 1.2 Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 14)
plot_pc_variance_explained(cds)

## check batch
# plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE,reduction_method = "PCA")
plot_cells(cds, color_cells_by="subtype", 
           label_cell_groups=FALSE,reduction_method = "PCA")
# plot_cells(cds, color_cells_by="data.sets", label_cell_groups=FALSE,reduction_method = "PCA")

#### Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "datasetsID", num_dim = 100)
# plot_cells(cds, color_cells_by="datasetsID", label_cell_groups=FALSE,reduction_method = "Aligned")

#### Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, preprocess_method = "Aligned")
# p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="subtype") + ggtitle('cds.umap')
# p1

#### Step 4: umap from Seurat data
cds.embed <- cds@int_colData$reducedDims$PCA
int.embed <- Embeddings(NE.data, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="subtype") +
  plot_cells(cds, reduction_method="UMAP", color_cells_by="datasetsID") 
p2

#### Step 5: Monocle3 Cluster
cds <- cluster_cells(cds, resolution=0.005)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2);p

#### Step 6: Identification trajectory
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "subtype",
           label_groups_by_cluster=T,
           label_leaves=T,
           label_branch_points=T)

#### Step 7: Order
# check# order_cells
FeaturePlot(NE.data,features = NEP100,
            reduction = "umap",pt.size = 0.1)&NoAxes()

cds <- order_cells(cds) 

col1 <- c('#b30c2a','#ce5c69','#e0a980','#f4c889','#a281b1','#5770a6','#364888')

rbPal <- colorRampPalette(col1)


plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE,  label_branch_points = F,label_roots = F,
           show_trajectory_graph = F)&NoAxes()&ggtitle("")&
  ggplot2::scale_colour_gradientn(name = "Pseudotime", 
                                  colours = rev(rbPal(100)))

#########################################
##################Figure 6D##############
#########################################
library(ggplot2)
library(ggrepel)
library(ggnewscale)
data <- read.csv('output_memento.csv')
# 根据FDR值设置显著性标记
data$significance <- ifelse(data$de_pval < 0.01, "de_pval < 0.01", "de_pval ≥ 0.01")
# 将cluster转换为因子类型
data$subtype <- factor(data$subtype,levels = c('NEC_O','Prol_H','Prol_E','EMT_R'))
data <- subset(data,data$de_coef>0)
# 计算每个subtype的最小和最大logFC值
col <- data %>%
  group_by(subtype) %>%
  summarise(min_logFC = min(de_coef), max_logFC = max(de_coef))
# 获取每个subtype中绝对值最大的前10个logFC
mark <- subset(data,data$gene %in% feature)
mark <- data %>%  group_by(subtype) %>%  top_n(10, abs(de_coef))
# 导入所需的库
# 创建 ggplot 图形
p <- ggplot() +
  geom_col(data = col,
           aes(x = subtype, y = min_logFC),
           fill = "grey90", alpha = 0.4, width = 0.8) +
  geom_col(data = col,
           aes(x = subtype, y = max_logFC),
           fill = "grey90", alpha = 0.4, width = 0.8) +
  geom_jitter(data = subset(data, abs(de_coef) > 0.5),
              aes(x = subtype, y = de_coef, group = subtype, color = significance),
              size = 0.8, width = 0.4) +
  geom_text_repel(data = mark,
                  aes(x = subtype, y = de_coef, label = gene),
                  size = 3.5, max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  color = 'black', force = 1) +
  geom_tile(data = mark,
            aes(x = subtype, y = 0, fill = subtype),
            height = 0.8, color = "black", alpha = 0.6, show.legend = FALSE) +
  scale_y_continuous(breaks = seq(-4, 4, 2), expand = expansion(add = c(0.1, 0.1))) +
  
  # 手动设置显著性标记的颜色
  scale_color_manual(values = c("de_pval ≥ 0.01" = "black",
                                "de_pval < 0.01" = "#ce5c69")) +

  scale_fill_manual(values = c("NEC_O" = "#5770a6",
                               "Prol_H" = "#a281b1",
                               "Prol_E" = "#f4c889",
                               "EMT_R" = "#ce5c69")
  ) +

  guides(color = guide_legend(override.aes = list(size = 3))) +
  ggnewscale::new_scale_color() +
  scale_colour_manual(values = c("NEC_O" = "white",
                                 "Prol_H" = "white",
                                 "Prol_E" = "white",
                                 "EMT_R" = "white")) +
  geom_text(data = mark,
            aes(x = subtype, y = 0, label = subtype, color = subtype),
            size = 5, show.legend = FALSE) +
  
  labs(x = "Subtype", y = "de_coef") +
  
  theme_minimal() +
  theme(axis.title = element_text(size = 20),
        axis.text.y = element_text(size = 16),
        axis.line.y = element_line(size = 1),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = c(0.11, 0.95),
        legend.title = element_blank(),
        legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 14),
        panel.grid = element_blank())
p

ggsave("multi-comparison dot plot.pdf", plot = p, width = 10, height = 6, device = cairo_pdf)

#########################################
##################Figure 7F-G##############
#########################################

data <- NE_sce@meta.data
table(data$sampletype)
table(data$sampleID)

count_data <- data %>%
  group_by(sampleID) %>%
  summarise(count = n())  # 计算每个 sampleID 的出现次数

# 筛选出计数大于等于 100 的 sampleID
valid_sampleIDs <- count_data %>%
  filter(count >= 50) %>%
  pull(sampleID)  # 提取符合条件的 sampleID

# 使用筛选出的 sampleID 过滤原始数据
data <- data %>%
  filter(sampleID %in% valid_sampleIDs)



# 安装并加载必要的包
# install.packages("pacman") # 如果没有pacman，请先安装
library(pacman)
p_load(openxlsx, ggplot2, dplyr)

data <- read.csv('data.csv')
# 设置因子水平
data <- data %>%
  filter(!is.na(subtype_total2)) %>%  # 过滤掉 subtype_total2 为 NA 的行
  group_by(sampleID, sampletype, subtype_total2) %>%  # 包含 sampletype 列
  summarise(count = n(), .groups = 'drop') %>%  # 计算每组的计数
  group_by(sampleID, sampletype) %>%  # 再次按 sampleID 和 sampletype 分组以计算比例
  mutate(percentage = (count / sum(count)) * 100) %>%  # 计算比例
  ungroup() %>%  # 解除分组
  arrange(sampleID, subtype_total2)  # 按 sampleID 和 subtype_total2 排序

data$sampletype <- factor(data$sampletype, levels = c("mCSPC", "NEPC"))
data$subtype_total2 <- factor(data$subtype_total2, levels = c('CNE','Prol1','Prol2','EMT'))


number <- max(data$id)
angle <- 90 - 360 * (data$id - 0.5) / number
data$hjust <- ifelse(angle < -90, 1, 0)  # 标签对齐方式
data$angle <- ifelse(angle < -90, angle + 180, angle)  # 翻转角度

# 为基线准备数据框
base_data <- data %>%
  group_by(sampletype) %>%
  summarize(start = min(id), end = max(id) - 3) %>%
  rowwise() %>%
  mutate(title = mean(c(start, end)))

base_data$angle <- c(-50, 90)

# 为网格（刻度）准备数据框
grid_data <- base_data
grid_data$end <- grid_data$end[c(nrow(grid_data), 1:(nrow(grid_data) - 1))] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# 创建径向柱状图
ggplot(data) +
  geom_bar(aes(x = as.factor(id), y = percentage, fill = subtype_total2), stat = "identity") +
  coord_radial(start = 0, end = 1 * pi) +
  labs(title = "", subtitle = "") +
  labs(x = "", y = "", fill = "") +
  scale_fill_manual(na.translate = FALSE,
                    values = c("CNE" = "#5770a6",
                               "Prol1" = "#a281b1",
                               "Prol2" = "#f4c889",
                               "EMT" = "#ce5c69")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(limits = c(-50, 110), breaks = c(0,  50,  100), labels = c(0,  50,  100)) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(color = "black", size = 18),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 2, "cm"),
    legend.text = element_text(size = 20),
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank(),
    panel.background = element_blank()
  ) +
  geom_text(data = data,
            aes(x = as.factor(id), y = 105, label = sampleID, hjust = hjust),
            size = 7, color = "grey20", angle = data$angle,
            inherit.aes = FALSE) +
  # 添加基线
  geom_segment(data = grid_data,
               aes(x = end, y = 100, xend = start, yend = 100),
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data,
               aes(x = end, y = 75, xend = start, yend = 75),
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data,
               aes(x = end, y = 50, xend = start, yend = 50),
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data,
               aes(x = end, y = 25, xend = start, yend = 25),
               colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +
  geom_segment(data = grid_data,
               aes(x = min(data$id) - 1, y = 0,
                   xend = min(data$id) - 1, yend = 100),
               alpha = 1, size = 0.6, inherit.aes = FALSE) +
  annotate("text",
           x = 0,
           y = 50,
           label = "Species\n threatened (%)",
           size = 6,
           angle = 90,
           hjust = 0.5) +
  #geom_segment(data = data,
  #             aes(x = id - 0.4, y = perBest.y * 100,
  #                 xend = id + 0.4, yend = perBest.y * 100),
  #             colour = "blue", size = 0.6, inherit.aes = FALSE) +
  geom_segment(data = base_data,
               aes(x = start, y = -5, xend = end, yend = -5),
               colour = "black", size = 1, inherit.aes = FALSE) +
  geom_text(data = base_data,
            aes(x = title, y = -15, label = sampletype),
            hjust = c(0.5, 0.5),
            colour = "black", size = 6, angle = base_data$angle,
            inherit.aes = FALSE)

p
#########################################
##################Figure 7H##############
#########################################
deg <- read.csv('output_memento.csv')

table(deg$subtype)

deg <- subset(deg,deg$de_coef>0)

mark <- deg %>%  group_by(subtype) %>%  top_n(50, abs(de_coef))

NEC_O_marker <- subset(mark$gene,mark$subtype=='NEC_O')
Prol_H_marker <- subset(mark$gene,mark$subtype=='Prol_H')
Prol_E_marker <- subset(mark$gene,mark$subtype=='Prol_E')
EMT_R_marker <- subset(mark$gene,mark$subtype=='EMT_R')


ne.data <- read.csv('G:/importance/undergraduated/数据集整理/bulk data（clinical）/NEPC/NE.csv',row.names = 1)

saveRDS(expr,'NE.expr.group.rds')

marker <- list(NEC_O=NEC_O_marker,
               Prol_H=Prol_H_marker,
               Prol_E=Prol_E_marker,
               EMT_R=EMT_R_marker)

library(IOBR)

expr <- read.csv("G:/importance/undergraduated/肾癌数据集/bulk/TCGA/TCGA.csv",row.names = 1)

rs2 <- IOBR::calculate_sig_score(pdata           = NULL,
                                 eset            = expr,
                                 signature       = list(intersection_up_sample),
                                 method          = "ssgsea",
                                 mini_gene_count = 0)

write.csv(rs2,'rs2.csv')

unlist(sig_group)

rs2 <- as.data.frame(rs2)
rownames(rs2) <- rs2[,1]
rs2 <- rs2[,-1]
rs2 <- t(rs2)
set.seed(123)


library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-0.4,0,0.4), c('#5770a6','white',"#ce5c69"))

set.seed(520)  # 为了可重复性
kmeans_result <- kmeans(t(rs2[1:4, ]), centers = 5)  # 注意转置，因列是样本

# 提取聚类分组信息
column_clusters <- kmeans_result$cluster
#column_clusters[column_clusters == 2] <- 4  # 将 2 替换为 4
#column_clusters[column_clusters == 5] <- 2  # 将 5 替换为 2（如果有的话）
column_clusters2 <- read.csv('bulk_group.csv')
column_clusters2$x <- factor(column_clusters2$x,levels = c('NEC-NE','Mixed-NE','EMT-NE','Low-spec'))


res <- as.data.frame(res)
row.names(group) <- group[,1]

res2 <- res[,-c(1,18)]
res2 <- as.matrix(res2)
Heatmap(res2,  # 不包含 NE 列
        name = "Z-score",
        #col = col_fun,
        na_col = "grey",
        border = TRUE,
        border_gp = gpar(lty = 1, lwd = 2, col = "black"),
        rect_gp = gpar(col = "white"),
        cluster_rows = F,#km = 4,
        cluster_columns = T,#column_km = 5,
        #row_names_side = "right",
        #row_title = "",
        row_split = group$subtype,
        # column_split = group$subtype,
        #row_title_side = "right",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_names = T,
        show_column_names = F,
        column_title = "",  # 添加列标题
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_gap = unit(10, "mm"),  # 调整行分割的宽度
        column_gap = unit(3, "mm"))  # 调整列分割的宽度
# right_annotation = ha_right,
#top_annotation = ha_top_combined)

#########三元图################
devtools::install_local("G:/迅雷下载/scTernary-main")  
data_exp_mat <- read.csv("G:/importance/undergraduated/数据集整理/bulk data（clinical）/NEPC/NE.csv",row.names = 1)
anno_signature_genes <- read.csv('bulk_markers.csv')
library(scTernary)
# 调用函数生成三元图数据
data_for_ternary <- generate_data_for_ternary(
  data_exp_mat = data_exp_mat,                  # 表达矩阵，包含基因表达数据
  anno_signature_genes = anno_signature_genes,  # 注释的标志基因列表
  gene_name_col = "Markers",                     # 基因名称列的名称
  gene_type_col = "subtype2",                  # 基因类型列的名称
  weight_by_gene_count = T,                   # 是否根据基因计数加权
  cutoff_exp =1,                               # 表达量的截止值
  prior_count = 1                      # 先验计数，用于平滑处理
)

cli <- read.csv('bulk_group.csv')
data_for_ternary <- as.data.frame(data_for_ternary)
data_for_ternary$ID <- rownames(data_for_ternary)

data_for_ternary <- left_join(data_for_ternary,cli,by='ID')


write.csv(data_for_ternary,'data_for_ternary.csv')
data_for_ternary <- read.csv('data_for_ternary.csv')
# 绘制三元图
vcdTernaryPlot(
  data_for_ternary,                            # 输入的数据，通常是经过预处理的三元图数据
  order_colnames = c(3, 1, 2),                # 列名的顺序，指定三元图的三个维度
  point_size = 1,                           # 点的大小
  group = data_for_ternary$subtype,  # 用于分组的变量，通常是 Seurat 聚类结果
  show_legend = TRUE,                         # 是否显示图例
  scale_legend = 0.8,                         # 图例的缩放比例
  legend_position = c(0.2, 0.5),              # 图例的位置，使用相对坐标（x, y）
  legend_vertical_space = 1,                   # 图例中各项之间的垂直间距
  legend_text_size = 1                         # 图例文本的大小
)


#########################################
##################Figure 7K##############
#########################################
library(limma)
library(oncoPredict)
library(parallel)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)
library(dplyr)
library(tibble)
library(ggpubr)
library(Hmisc)
library(tidyr)

CTRP2_Expr <- readRDS(file="G:/importance/undergraduated/12_R包开发/EZGen/inst/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res <- readRDS(file = "G:/importance/undergraduated/12_R包开发/EZGen/inst/DataFiles/Training Data/CTRP2_Res.rds")

#设置参数
batchCorrect<-"eb"
powerTransformPhenotype<-TRUE
removeLowVaryingGenes<-0.2
removeLowVaringGenesFrom<-'rawData'
minNumSamples=10
selection<- 1
printOutput=TRUE
pcr=FALSE
report_pc=FALSEcc=FALSE
rsq=FALSE
percent=80

expr <- read.csv("D:/下载/WeChat Files/wxid_4tmvpn3qli1m22/FileStorage/File/2025-05/geneMatrix.csv",row.names = 1)
expr <- as.matrix(expr)
#计算并自动输出
max-ppsize=500000
colnames(expr)
expr1 <- expr[,1:2]
expr2<- expr[,181:373]
options(expressions = 500000) 
Cstack_info()
memory.limit(size=8000000)
str.default()
object.size(expr1)
memory.profile()

calcPhenotype(trainingExprData=CTRP2_Expr,
              trainingPtype=CTRP2_Res,
              testExprData=expr1,
              batchCorrect=batchCorrect,
              powerTransformPhenotype=powerTransformPhenotype,
              removeLowVaryingGenes=removeLowVaryingGenes,
              minNumSamples=minNumSamples,
              selection=selection,
              printOutput=printOutput,
              pcr=pcr,
              removeLowVaringGenesFrom=removeLowVaringGenesFrom,
              report_pc=report_pc,
              percent=percent,
              rsq=rsq)

rs <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)

col_fun = colorRamp2(c(-2,0,2), c('#5770a6','white',"#ce5c69"))
rs <- t(rs)
rs <- Zscore(rs)

Zscore <- function(expr=av.exp,feature=feature){
  #expr2 <- log2(expr+1)
  #expr2<- expr2[which(row.names(expr2)%in% feature),]
  z_scores <-  t(scale(t(expr)))
  z_scores_normalized <- pmax(pmin(z_scores, 2), -2)
  z_scores_normalized <- as.matrix(z_scores_normalized)
  return(z_scores_normalized)
}

Heatmap(rs,  # 不包含 NE 列
        name = "Z-score",
        col = col_fun,
        na_col = "grey",
        border = TRUE,
        border_gp = gpar(lty = 1, lwd = 2, col = "black"),
        rect_gp = gpar(col = "white"),
        cluster_rows = T,#km = 4,
        cluster_columns = F,#column_km = 5,
        #row_names_side = "right",
        #row_title = "",
        #row_split = row_split,
        column_split = column_clusters2$x,
        #row_title_side = "right",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_names = T,
        show_column_names = F,
        column_title = "",  # 添加列标题
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_gap = unit(10, "mm"),  # 调整行分割的宽度
        column_gap = unit(3, "mm"))  # 调整列分割的宽度
# right_annotation = ha_right,
#top_annotation = ha_top_combined)
colnames(column_clusters2) <- c('Sample','Cluster')
rs <- as.data.frame(rs)
rs <- t(rs)
rs_transposed <- rownames_to_column(rs, var = "Sample")
rs_long <- rs_transposed %>%
  pivot_longer(-Sample, names_to = "Drug", values_to = "Value") %>%
  left_join(column_clusters2, by = "Sample")

# 计算每个分组的正数比例和平均值
summary_data <- rs_long %>%
  group_by(Cluster, Drug) %>%  # 按分组和药物分组
  summarise(
    PositiveRatio = mean(Value > 0),  # 计算正数比例
    MeanValue = mean(Value, na.rm = TRUE)  # 计算平均值
  ) %>%
  ungroup()  # 取消分组，以便后续操作

sedrug <- c('carboplatin','oxaliplatin',#Platinum
            'etoposide','teniposide','topotecan','doxorubicin',#TOP
            'bendamustine','bleomycin.A2','chlorambucil','dacarbazine','gemcitabine',#DNA damage
            'clofarabine','cytarabine.hydrochloride','fluorouracil','methotrexate',#Antimetabolite
            'navitoclax','ABT.199',#BCL2
            'alisertib','barasertib',#AURK
            'olaparib','veliparib')#PARP
summary_data2 <- subset(summary_data,summary_data$Drug%in% sedrug)
rs2 <- rs[,colnames(rs)%in% sedrug]


summary_data2$Drug <- factor(summary_data2$Drug,levels = sedrug)

write.csv(rs2,'drug.csv')
write.csv(summary_data2,'drug_summary.csv')
# 绘制气泡图
ggplot(subset(summary_data2,!summary_data2$Cluster=='Low-spec'), aes(x = Cluster, y = Drug, size = PositiveRatio, color = MeanValue)) +
  geom_point(alpha = 1) +  # 设置气泡透明度
  scale_size_continuous(range = c(2, 8)) +  # 设置气泡大小范围
  scale_color_gradientn(colors = c('#5770a6', 'white', '#ce5c69')) +  # 设置颜色渐变
  theme_minimal() +  # 使用简约主题
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
        panel.spacing = unit(2, "lines"), strip.background = element_blank())+
  
  labs(title = "",
       x = "Subtype",
       y = "Drug",
       size = "Resistance Ratio",
       color = "IC50(Z-score)")


