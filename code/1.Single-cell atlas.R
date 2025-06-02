####################single-cell QC###################################
library(DoubletFinder)
library(tidyverse)
library(Seurat)
library(patchwork)
library(dplyr)
library(future)
plan("multisession", workers = 2)
plan()
options(future.globals.maxSize = 100 * 1024^3)
rm(list = ls());gc()#删除所有环境变量整理内存空间


#文件夹格式
root_directory <- "G:/importance/undergraduated/数据集整理/single cell/GSE250189（验证）/"  # 对于 Windows 系统
folders <- list.dirs(root_directory, full.names = TRUE, recursive = FALSE)

#csv格式
folder_path <- "G:/importance/undergraduated/数据集整理/single cell/GSE210358"
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = F)

#txt格式
folder_path <- "G:/importance/undergraduated/数据集整理/single cell/GSE137829"
file_list <- list.files(path = folder_path, pattern = "\\.txt$", full.names = F)

for (i in seq_along(folders)) {
  #matrix_data <- read.csv(file.path(folder_path, file_list[i]))  # 使用完整路径读取文件
  #matrix_data <- t(matrix_data)
  #allcell.data <- as.data.frame(matrix_data)
  #colnames(allcell.data) <- as.character(allcell.data[1, ])
  #allcell.data <- allcell.data[-1,]
  
  #expr <- matrix_data %>%
  #  group_by(across(1)) %>%  # 根据第一列进行分组
  #  summarize(across(everything(), max, .names = "{col}"))  # 计算每列的最大值

  allcell.data <- Read10X(data.dir =folders[i])
  #GSE240056 <- read.table('G:/importance/undergraduated/数据集整理/single cell/GSE240056（验证）/GSE240056counts.txt', sep="\t", header=T)  # 使用完整路径读取文件
  allcell <- CreateSeuratObject(counts = allcell.data, project = "GSE250189",min.cells = 0, min.features = 0)

  sampleID <- unlist(strsplit(basename(folders[i]), "_"))[1]

  
  allcell@meta.data$sampleID <- sampleID
  
  # 标准化表达数据与前期准备工作
  allcell <- NormalizeData(allcell, normalization.method = "LogNormalize", scale.factor = 10000)
  allcell <- FindVariableFeatures(allcell, selection.method = "vst", nfeatures = 2000)
  allcell <- ScaleData(allcell, features = rownames(allcell))
  allcell <- RunPCA(allcell, features = VariableFeatures(object = allcell))
  allcell <- RunUMAP(allcell, dims = 1:20)
  allcell <- FindNeighbors(allcell, dims = 1:20) %>% FindClusters(resolution = 0.6)
  
  # 预测双联细胞
  sweep.res.list <- paramSweep(allcell, PCs = 1:20)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
  
  DoubletRate <- ncol(allcell) * 8 * 1e-6
  homotypic.prop <- modelHomotypic(allcell$seurat_clusters)
  
  nExp_poi <- round(DoubletRate * nrow(allcell@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  allcell <- doubletFinder(allcell, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
  
  # 计算线粒体、核糖体和热休克蛋白的表达比例
  mt.genes <- rownames(allcell)[grep("^(MT-|MTRNR|MT\\.)", rownames(allcell))]
  C <- GetAssayData(object = allcell, slot = "counts")
  percent.mt <- colSums(C[mt.genes,]) / Matrix::colSums(C) * 100
  allcell <- AddMetaData(allcell, percent.mt, col.name = "percent.mt")
  
  rb.genes <- rownames(allcell)[grep("^RP[SL]", rownames(allcell))]
  percent.rb <- colSums(C[rb.genes,]) / Matrix::colSums(C) * 100
  allcell <- AddMetaData(allcell, percent.rb, col.name = "percent.rb")
  
  hsp.genes <- rownames(allcell)[grep("^(HSP|DNAJ())", rownames(allcell))]
  percent.hsp <- colSums(C[hsp.genes,]) / Matrix::colSums(C) * 100
  allcell <- AddMetaData(allcell, percent.hsp, col.name = "percent.hsp")
  
  hb.genes <- rownames(allcell)[grep("^HB[^(P)]", rownames(allcell))]
  percent.hb <-  colSums(C[hb.genes,]) / Matrix::colSums(C) * 100
  allcell <- AddMetaData(allcell, percent.hb, col.name = "percent.hb")
  
  
  VlnPlot(allcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.rb', 'percent.hsp','percent.hb'), ncol = 3,group.by = "orig.ident")
  ggsave(filename = paste0('./QC VlnPlot/',sampleID,' QC VlnPlot.pdf'),width = 12,height = 8)
  
  # 去除预测双联细胞、低质量细胞(nFeature_RNA > 250, nCount_RNA > 500, percent.mt < 20)与红细胞(percent.hb < 5)
  allcell.qc1 <- subset(allcell, subset = nFeature_RNA > 250 & nCount_RNA > 500 & percent.mt < 20 & percent.hb < 5 &!!sym(names(allcell@meta.data[8])) == "Singlet")

  # 去除低表达基因(<3 cells)
  count.data <- GetAssayData(object = allcell.qc1, slot = "counts")
  
  nonzero <- count.data > 0
  
  low.genes <- Matrix::rowSums(nonzero) <= 3
  
  allcell.qc2 <- allcell.qc1[!low.genes, ]
  
  # 避免意外的噪音和解离的表达假象，排除了与线粒体（50个基因）热休克蛋白（178个基因）核糖体（1,253个基因）相关的基因

  genes_to_remove <- c(mt.genes, hsp.genes, rb.genes )
  
  allcell.qc3 <- allcell.qc2[!rownames(allcell.qc2) %in% genes_to_remove, ]

  allcell.qc3@assays$RNA$data <-  NULL
  allcell.qc3@assays$RNA$scale.data <-  NULL
  
  saveRDS(allcell.qc3, paste0(root_directory,sampleID,'.rds'))
}

####################annotation##################

folder_path <- "G:/importance/undergraduated/数据集整理/single cell/暂存/"
file_list.path <- list.files(path = folder_path, pattern = "\\.rds$", full.names = T)

sampletype <- read.csv('G:/importance/undergraduated/数据集整理/single cell/暂存/sampletype.csv')

for (i in 1:length(file_list.path)) {
  print(i)
  seurat.data <- read_rds(file_list.path[i])
  seurat.data@assays$RNA$data <-  NULL
  seurat.data@assays$RNA$scale.data <-  NULL
  names(seurat.data@meta.data)[1] <-'datasetsID'
  seurat.data@meta.data$sampletype <- sampletype$sampletype[match(seurat.data@meta.data$sampleID, sampletype$sampleID)]
  saveRDS(seurat.data, file_list.path[i])
  rm(seurat.data);gc()
}

for (i in 1:length(file_list.path)) {
  seurat.data <- read_rds(file_list.path[1])
  seurat.data@meta.data[["sampletype"]] <-'Primary'
  saveRDS(seurat.data, file_list.path[1])
  rm(seurat.data);gc()
}

table(allcell@meta.data[["orig.ident"]])

allcell@meta.data$sampleID <- sub("-.*", "", allcell@meta.data[["orig.ident"]])
allcell@meta.data$sampletype <- ifelse(grepl("TP3", allcell@meta.data[["orig.ident"]]), "early",
                                       ifelse(grepl("TP4|TP5", allcell@meta.data[["orig.ident"]]), "transition",
                                              ifelse(grepl("TP6", allcell@meta.data[["orig.ident"]]), "endpoint", NA)))


for (i in 1:16) {
  count.data <- GetAssayData(object =GSE215943 , slot = "counts")
  seurat.data@meta.data[["sampletype"]] <-'Primary'
  saveRDS(seurat.data, file_list.path[1])
  
}

library("biomaRt")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "useast")
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

GSE215943@meta.data <- GSE215943@meta.data %>% 
  mutate(celltype = case_when(
    sampleID == 'GSM6647951' ~ 'DMSO',
    sampleID == 'GSM6647952' ~ 'ENZ_D4',
    sampleID == 'GSM6647953' ~ 'ENZ_D7',
    sampleID == 'GSM6647954' ~ 'ENZ_D14'
  ))
##################Combination of multiple samples################
# 读入和处理每个样本的数据
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
library(harmony)
library(ROGUE)
library(DoubletFinder)
library(DESeq2)
library(tibble)

BiocManager::install("SingleR")

rm(list = ls());gc()
nbrOfWorkers()
plan("multisession", workers = 1)
options(future.globals.maxSize = 40 * 1024^3)

folder_path <- "G:/importance/undergraduated/数据集整理/single cell/allsample"
file_list.path <- list.files(path = folder_path, pattern = "\\.rds$", full.names = T)
file_list.name <- list.files(path = folder_path, pattern = "\\.rds$", full.names = F)
file_ids <- sub("\\.rds$", "", file_list.name)

seurat.data.list <- list()
for (i in 1:length(file_list.path)) {
  print(i) 
  seurat.data.list[i] <- read_rds(file_list.path[i])
}

#Combination
merge.seurat.data <- merge(seurat.data.list[[1]],
                           y=c(seurat.data.list[[2]], seurat.data.list[[3]], seurat.data.list[[4]], seurat.data.list[[5]], seurat.data.list[[6]], seurat.data.list[[7]], seurat.data.list[[8]], seurat.data.list[[9]], seurat.data.list[[10]], seurat.data.list[[11]], seurat.data.list[[12]], seurat.data.list[[13]], seurat.data.list[[14]], seurat.data.list[[15]], seurat.data.list[[16]], seurat.data.list[[17]], seurat.data.list[[18]], seurat.data.list[[19]], seurat.data.list[[20]], seurat.data.list[[21]], seurat.data.list[[22]], seurat.data.list[[23]], seurat.data.list[[24]], seurat.data.list[[25]], seurat.data.list[[26]], seurat.data.list[[27]], seurat.data.list[[28]], seurat.data.list[[29]], seurat.data.list[[30]], seurat.data.list[[31]], seurat.data.list[[32]], seurat.data.list[[33]], seurat.data.list[[34]], seurat.data.list[[35]], seurat.data.list[[36]], seurat.data.list[[37]], seurat.data.list[[38]], seurat.data.list[[39]], seurat.data.list[[40]], seurat.data.list[[41]], seurat.data.list[[42]], seurat.data.list[[43]], seurat.data.list[[44]], seurat.data.list[[45]], seurat.data.list[[46]], seurat.data.list[[47]], seurat.data.list[[48]], seurat.data.list[[49]], seurat.data.list[[50]], seurat.data.list[[51]], seurat.data.list[[52]], seurat.data.list[[53]], seurat.data.list[[54]], seurat.data.list[[55]], seurat.data.list[[56]], seurat.data.list[[57]], seurat.data.list[[58]], seurat.data.list[[59]], seurat.data.list[[60]], seurat.data.list[[61]], seurat.data.list[[62]], seurat.data.list[[63]], seurat.data.list[[64]], seurat.data.list[[65]], seurat.data.list[[66]], seurat.data.list[[67]], seurat.data.list[[68]], seurat.data.list[[69]], seurat.data.list[[70]]),
)
#Eliminate batch effect

RemoveBatch = function(seurat.data, 
                       batchID = "Data.sets", 
                       methods = "harmony", n.pcs = 50){
  seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
    ScaleData(verbose = F) %>%
    RunPCA(npcs = n.pcs, verbose = F)
  
  seurat.data <- seurat.data %>% RunHarmony(batchID, plot_convergence = T)
  
  seurat_int <- seurat.data %>% 
    RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = F)
  
  return(seurat_int)
}

seurat.harmony <- RemoveBatch(seurat.data = seurat.harmony, n.pcs = 50,methods = "harmony",batchID = "datasetsID" )

#Clustering and grouping
#selection parameter
ElbowPlot(seurat.harmony, ndims=50, reduction="nCount_RNA") 
pct <- seurat.harmony[["pca"]]@stdev / sum(seurat.harmony[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5),sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
#pc.use <- cumu > 90
ElbowPlot(seurat.harmony)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

seurat.harmony <- FindNeighbors(seurat.harmony,reduction = 'pca',dims = 1:pc.use)#选择多少个PCA维度进行分析

seurat.harmony <- FindClusters(seurat.harmony,resolution = seq(from = 0.1,to = 0.5, by = 0.1))#测试选择分辨率大小

clustree(seurat.harmony) 

seurat.harmony <- FindClusters(seurat.harmony,resolution = 0.3)

#UMAP
seurat.harmony <- RunUMAP(seurat.harmony, 
                          reduction = "harmony", 
                          dims = 1:pc.use, 
                          min.dist = 0.5, 
                          n.neighbors = 300, 
                          seed.use = 42)

DimPlot(seurat.harmony, reduction = "umap", group.by = "datasetsID", pt.size=1.4)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank())


#检验细胞簇纯度，引用张泽民教授
seurat.harmony <- JoinLayers(seurat.harmony)
sce_sub <- seurat.harmony %>% subset(., downsample = 100)
expr <- GetAssayData(sce_sub, slot = 'counts') %>% as.matrix()
meta <- sce_sub@meta.data

rogue(expr = expr, 
      labels = meta$seurat_clusters, 
      samples = meta$sampleID, 
      platform = "UMI",
      span = 0.9) %>% rogue.boxplot()

#手动分群#
features = c(
  # Epithelial
  'KRT8', 'CDH1', 'DST','LCN2','EPCAM',
  # Luminal
  'AR', 'KLK2', 'KLK3', 'TMPRSS2', 
  # NE
  'SYP', 'CHGA', 'CHGB', 'ENO2',
  # Fib
  'LUM', 'DCN', 'FBLN1','COL1A1'  ,'COL6A3',
  # MyoFib
  'ACTA2', 'RGS5', 'NOTCH3', 'PDGFRB', 'NDUFA4L2',
  # Ecs
  'CLDN5', 'PECAM1', 'VWF',  'CDH5','CD200',
  # Mast
  'PTPRC','KIT', 'TPSB2', 'MS4A2', 'TPSAB1',
  # Mye
  'CD68', 'HLA-DQB1', 'APOE', 'C1QA',
  # T
  'CD3D', 'CD3E', 'IL7R', 'NKG7',
  # NK
  'GZMB', 'KLRD1',
  # B
  'CD19', 'MS4A1', 'CCR7', 'CD79A',
  # Plasma
  'JCHAIN', 'MZB1', 'ZBP1', 'IGJ'
  
)

DotPlot(object = seurat.harmony,features =features ,scale.by = "size") + coord_flip()
#手动注释
Epithelial_cluster <- c()
Luminal_cluster <- c()
NE_cluster<- c()
Ecs_cluster <-c()
Fib_cluster <- c()
MyoFib_cluster<- c()
Mast_cluster <- 
Mye_cluster <-c() 
T_cluster<- c()
NK_cluster <- 
B_cluster<- c()
Plasma_cluster <- 

seurat.harmony@meta.data <- seurat.harmony@meta.data %>% 
  mutate(celltype = case_when(
    seurat_clusters %in% Plasma_cluster ~ 'Plasma',
    seurat_clusters %in% B_cluster ~ 'B',
    seurat_clusters %in% NK_cluster ~ 'NK',
    seurat_clusters %in% T_cluster ~ 'T',
    seurat_clusters %in% Mye_cluster ~ 'Mye',
    seurat_clusters %in% Fib_cluster ~ 'Fib',
    seurat_clusters %in% MyoFib_cluster ~ 'MyoFib',
    seurat_clusters %in% Mast_cluster ~ 'Mast',
    seurat_clusters %in% Ecs_cluster ~ 'Ecs',
    seurat_clusters %in% Epithelial_cluster ~ 'Epithelial',
    seurat_clusters %in% Luminal_cluster ~ 'Luminal',
    seurat_clusters %in% NE_cluster ~ 'NE'))

seurat.harmony@meta.data$celltype <- factor(seurat.harmony@meta.data$celltype,
                                            levels =  c('Epithelial','Luminal', 'NE',
                                                                'Fib','MyoFib','Ecs','Mast','Mye',
                                                                'T','NK','B','Plasma'))
table(Idents(seurat.harmony) )
table(seurat.harmony@meta.data$seurat_clusters )
table(seurat.harmony@meta.data$celltype )
Idents(seurat.harmony) <- "celltype"


DimPlot(seurat.harmony, reduction = "umap", group.by = "celltype", pt.size=1.4)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank())


####################visualization##################
library(plot1cell)
library(Seurat) 
library(tidyverse)
library(stringr)
library(RColorBrewer)

# Set the color palette
cluster_colors<-c('#b30c2a','#ce5c69','#F99999','#e0a980','#f4c889','#bdd5a3','#86a979','#519981','#8ba1c6','#5770a6','#a281b1','#735c88')
datasetsID_colors<-c('#b30c2a','#ce5c69','#e0a980','#f4c889','#bdd5a3','#519981','#8ba1c6','#5770a6','#a281b1')
sampletype_colors <- c('#ce5c69','#f4c889','#bdd5a3','#5770a6','#a281b1')
NE_colors<-c('#ce5c69','#5770a6')
tissuesource_colors <- c('#b30c2a','#ce5c69','#e0a980','#f4c889','#bdd5a3','#519981','#8ba1c6','#5770a6','#a281b1')

circ_data <- prepare_circlize_data(seurat.harmony, scale = 0.8 )
set.seed(20000709)

circ_data_sorted <- circ_data[order(circ_data$celltype), ]
circ_data_sorted$Cluster <- factor(circ_data_sorted$Cluster,
                                   levels = c('Epithelial','Luminal', 'NE',
                                                                'Fib','MyoFib','Ecs','Mast','Mye',
                                                                'T','NK','B','Plasma'))
#UMAP plotting function
plot_circlize_change <- function (data_plot, do.label = T, contour.levels = c(0.2, 0.3), 
                                  pt.size = 0.5, kde2d.n = 1000, contour.nlevels = 100, bg.color = "#F9F2E4", 
                                  col.use = NULL, label.cex = 0.5, labels.cex = 0.5, circos.cex = 0.5 ,repel = FALSE) 
{
  centers <- data_plot %>% dplyr::group_by(Cluster) %>% summarise(x = median(x = x), 
                                                                  y = median(x = y))
  z <- MASS::kde2d(data_plot$x, data_plot$y, n = kde2d.n)
  celltypes <- names(table(data_plot$Cluster))
  cell_colors <- (scales::hue_pal())(length(celltypes))
  if (!is.null(col.use)) {
    cell_colors = col.use
    col_df <- data.frame(Cluster = celltypes, color2 = col.use)
    cells_order <- rownames(data_plot)
    data_plot <- merge(data_plot, col_df, by = "Cluster")
    rownames(data_plot) <- data_plot$cells
    data_plot <- data_plot[cells_order, ]
    data_plot$Colors <- data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0.01, 
                                                            0), track.height = 0.01, gap.degree = c(rep(2, (length(celltypes) - 
                                                                                                              1)), 12), points.overflow.warning = FALSE)
  circos.initialize(sectors = data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y = data_plot$dim2, 
               bg.border = NA, panel.fun = function(x, y) {
                 circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + 
                               mm_y(4), CELL_META$sector.index, cex = labels.cex,
                             col = "black", facing = "bending.inside", niceFacing = T)
                 #circos.axis(labels.cex = 0.3, col = "black", labels.col = "black")
                 circos.axis(labels.cex = circos.cex, col = "black", labels.col = "black")
               })
  for (i in 1:length(celltypes)) {
    dd <- data_plot[data_plot$Cluster == celltypes[i], ]
    circos.segments(x0 = min(dd$x_polar2), y0 = 10, x1 = max(dd$x_polar2), 
                    y1 = 10, col = cell_colors[i], lwd = 5,sector.index = celltypes[i])
  }
  text(x = 1, y = 0.1, labels = "Cluster", cex = 0.4, col = "black", 
       srt = -90)
  points(data_plot$x, data_plot$y, pch = 19, col = alpha(data_plot$Colors, 
                                                         0.2), cex = pt.size)
  contour(z, drawlabels = F, nlevels = 50, lty = 2,lwd = 2, levels = contour.levels, 
          col = "black", add = TRUE)
  if (do.label) {
    if (repel) {
      textplot(x = centers$x, y = centers$y, words = centers$Cluster, 
               cex = label.cex, new = F, show.lines = F)
    }
    else {
      text(centers$x, centers$y, labels = centers$Cluster, 
           cex = label.cex, col = "black")
    }
  }
}
#########################################
##################Figure 2A##############
#########################################
plot_circlize_change(circ_data_sorted,do.label =F, pt.size = 0.1, 
                     col.use = cluster_colors ,
                     bg.color = 'white', 
                     kde2d.n = 1000, 
                     repel = F, 
                     labels.cex = 1, 
                     circos.cex = 0.5,
                     label.cex = 1)
add_track(circ_data_sorted, 
          group = "datasetsID", track_lwd = 3,
          colors = datasetsID_colors, track_num = 2) 
add_track(circ_data_sorted, 
          group = "sampletype",track_lwd = 3,
          colors = sampletype_colors, track_num = 3)
add_track(circ_data_sorted, 
          group = "NE",track_lwd = 3,
          colors = NE_colors, track_num = 4)
add_track(circ_data_sorted, 
          group = "tissuesource",track_lwd = 3,
          colors = tissuesource_colors, track_num = 5)
###################datasetsID
umap_plot1 <- DimPlot(seurat.harmony, 
                      cols = cluster_colors, 
                      reduction = "umap", 
                      group.by = "celltype12", 
                      pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  guides(color = guide_legend(override.aes = list(shape = 15,size=10)))
legend <- get_legend(umap_plot1)
plot(legend)
DimPlot(seurat.harmony, 
        cols = datasetsID_colors, 
        reduction = "umap", 
        group.by = "datasetsID", 
        pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",plot.title = element_blank()) 
###################sampletype
umap_plot1 <- DimPlot(seurat.harmony, 
                      cols = sampletype_colors, 
                      reduction = "umap", 
                      group.by = "sampletype", 
                      pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
legend <- get_legend(umap_plot1)
plot(legend)
seurat.harmony@meta.data$sampletype <- factor(seurat.harmony@meta.data$sampletype,
                                              level =c('Primary', 'mCSPC', 'CRPC', 'mCRPC', 'NEPC'))
DimPlot(seurat.harmony, 
        cols = sampletype_colors, 
        reduction = "umap", 
        group.by = "sampletype", 
        pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",plot.title = element_blank()) 
###################NE
umap_plot1 <- DimPlot(seurat.harmony, 
                      cols = NE_colors, 
                      reduction = "umap", 
                      group.by = "NE", 
                      pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
legend <- get_legend(umap_plot1)
plot(legend)
DimPlot(seurat.harmony, 
        cols = NE_colors, 
        reduction = "umap", 
        group.by = "NE", 
        pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",plot.title = element_blank()) 
###################tissuesource
umap_plot1 <- DimPlot(seurat.harmony, 
                      cols = tissuesource_colors, 
                      reduction = "umap", 
                      group.by = "tissuesource", 
                      pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
legend <- get_legend(umap_plot1)
plot(legend)
DimPlot(seurat.harmony, 
        cols = tissuesource_colors, 
        reduction = "umap", 
        group.by = "tissuesource", 
        pt.size = 1) + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",plot.title = element_blank()) 

#########################################
##################Figure 2B##############
#########################################
features = c( 
  # Epithelial
  'KRT8', 'CDH1', 'DST','LCN2','EPCAM',
  # Luminal
  'AR', 'KLK2', 'KLK3', 'TMPRSS2', 
  # NE
  'SYP', 'CHGA', 'CHGB', 'ENO2',
  # Fib
  'LUM', 'DCN', 'FBLN1','COL1A1'  ,'COL6A3',
  # MyoFib
  'ACTA2', 'RGS5', 'NOTCH3', 'PDGFRB', 'NDUFA4L2',
  # Ecs
  'CLDN5', 'PECAM1', 'VWF',  'CDH5','CD200',
  # Mast
  'PTPRC','KIT', 'TPSB2', 'MS4A2', 'TPSAB1',
  # Mye
  'CD68', 'HLA-DQB1', 'APOE', 'C1QA',
  # T
  'CD3D', 'CD3E', 'IL7R', 'NKG7',
  # NK
  'GZMB', 'KLRD1',
  # B
  'CD19', 'MS4A1', 'CCR7', 'CD79A',
  # Plasma
  'JCHAIN', 'MZB1', 'ZBP1', 'IGJ'
  
)

DotPlot(seurat.harmony,features = features)+ coord_flip()
library(paletteer)

cluster_colors <- c('#b30c2a','#ce5c69','#F99999','#e0a980','#f4c889',
                    '#bdd5a3','#86a979','#519981','#8ba1c6','#5770a6',
                    '#a281b1','#735c88')
reversed_colors <- rev(cluster_colors)
repeats <- c(5, 4, 4, 5, 5, 5, 5, 4, 4, 2, 4, 4)
reversed_features <- rev(features)
color_vector <- rep(reversed_colors[1:length(repeats)], repeats)
#stack=T绘制堆叠小提琴图
VlnPlot(seurat.harmony,features=reversed_features,
        group.by="celltype",
        stack=T,cols=color_vector,flip=T,pt.size = 0
)+NoLegend()


#####################Calculate the score of the gene set：AddModuleScore###########
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

#####################celltype
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
###################sampletype
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
#########################################
##################Figure 2C##############
#########################################
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
#########################################
##################Figure 2D##############
#########################################
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
  #############scatter diagram
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
#########################################
##################Figure 2D##############
#########################################
for (i in 1:27) {
  datasetID <- colnames(data)[i]
  print(datasetID)
  
  comparisons_list <- list(
    c("NEPC", "CRPC"),
    c("NEPC", "mCRPC"),
    c("NEPC", "mCSPC") ,
    c("NEPC", "Primary") # 添加更多的比较组
  )
  
  ######Box plot
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


#########################################
##################Figure 2E##############
#########################################
# filtered_data 是一个数据框，包含基因名称
genes_to_plot <- filtered_data$gene

# 计算表达百分比
expression_data <- seurat.harmony@assays$RNA$data[genes_to_plot, ]  # 获取基因表达数据
expression_df <- as.data.frame(t(expression_data))  # 转置数据框

# 添加细胞类型信息
expression_df$celltype <- seurat.harmony$celltype12

# 计算每个基因在每个细胞类型中的表达百分比
percent_expression <- expression_df %>%
  group_by(celltype) %>%
  summarise(across(all_of(genes_to_plot), ~ mean(. > 0) * 100, .names = "percent_{col}")) %>%
  pivot_longer(-celltype, names_to = "gene", values_to = "percent")  # 转换为长格式

# 按表达百分比排序基因
percent_expression <- percent_expression %>%
  arrange(desc(percent))

# 获取排序后的基因列表
sorted_genes <- unique(percent_expression$gene)

sorted_genes <- sub('^\\percent_','',sorted_genes)

###Dotplot function
custom_dotplot <- function(
    object, 
    features, 
    assay = NULL, 
    cols = c("lightgrey", "blue"), 
    col.min = -2.5, 
    col.max = 2.5, 
    dot.min = 0, 
    dot.scale = 6, 
    idents = NULL, 
    group.by = NULL, 
    split.by = NULL, 
    cluster.idents = FALSE, 
    scale = TRUE, 
    scale.by = "radius", 
    scale.min = NA, 
    scale.max = NA) {
  # 设置默认assay
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  
  # 处理颜色和缩放函数
  split.colors <- !is.null(split.by) && !any(cols %in% rownames(brewer.pal.info))
  scale.func <- switch(scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
  # 处理features
  feature.groups <- NULL
  if (is.list(features) || any(!is.na(names(features)))) {
    feature.groups <- unlist(sapply(1:length(features), function(x) {
      rep(names(features)[x], each = length(features[[x]]))
    }))
    
    if (any(is.na(feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, immediate. = TRUE)
    }
    
    features <- unlist(features)
    names(feature.groups) <- features
  }
  
  # 获取细胞数据
  cells <- unlist(CellsByIdentities(object = object, cells = colnames(object[[assay]]), idents = idents))
  data.features <- FetchData(object = object, vars = features, cells = cells)
  
  # 处理分组信息
  data.features$id <- if (is.null(group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  
  # 确保id是因子类型
  if (!is.factor(data.features$id)) {
    data.features$id <- factor(data.features$id)
  }
  
  id.levels <- levels(data.features$id)
  data.features$id <- as.vector(data.features$id)
  
  # 处理split信息
  if (!is.null(split.by)) {
    splits <- FetchData(object = object, vars = split.by)[cells, split.by]
    
    if (split.colors) {
      if (length(unique(splits)) > length(cols)) {
        stop(paste0("Need to specify at least ", length(unique(splits)), " colors using the cols parameter"))
      }
      
      cols <- cols[1:length(unique(splits))]
      names(cols) <- unique(splits)
    }
    
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(splits)
    id.levels <- paste0(rep(id.levels, each = length(unique.splits)), "_", rep(unique.splits, times = length(id.levels)))
  }
  
  # 计算平均表达和表达百分比
  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, 1:(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, PercentAbove, threshold = 0)
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  
  names(data.plot) <- unique(data.features$id)
  
  # 处理聚类和排序
  if (cluster.idents) {
    mat <- do.call(rbind, lapply(data.plot, unlist))
    mat <- scale(mat)
    id.levels <- id.levels[hclust(dist(mat))$order]
  }
  
  # 整合数据用于绘图
  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$id <- x
    return(data.use)
  })
  
  data.plot <- do.call(rbind, data.plot)
  
  if (!is.null(id.levels)) {
    data.plot$id <- factor(data.plot$id, levels = id.levels)
  }
  
  # 检查分组数量
  ngroup <- length(levels(data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will not be scaled", call. = FALSE, immediate. = TRUE)
  } else if (ngroup < 5 && scale) {
    warning("Scaling data with a low number of groups may produce misleading results", call. = FALSE, immediate. = TRUE)
  }
  
  # 缩放平均表达
  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    
    if (scale) {
      data.use <- scale(log1p(data.use))
      data.use <- MinMax(data.use, min = col.min, max = col.max)
    } else {
      data.use <- log1p(data.use)
    }
    
    return(data.use)
  })
  
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  
  if (split.colors) {
    avg.exp.scaled <- as.numeric(cut(avg.exp.scaled, breaks = 20))
  }
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(data.plot$features.plot, levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  # 处理颜色
  if (split.colors) {
    splits.use <- unlist(lapply(data.plot$id, function(x) {
      sub(paste0(".*_(", paste(sort(unique(splits), decreasing = TRUE), collapse = "|"), ")$"), "\\1", x)
    }))
    
    data.plot$colors <- mapply(function(color, value) {
      colorRampPalette(colors = c("grey", color))(20)[value]
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  
  color.by <- if (split.colors) "colors" else "avg.exp.scaled"
  
  if (!is.na(scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  
  if (!is.na(scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  
  if (!is.null(feature.groups)) {
    data.plot$feature.groups <- factor(feature.groups[data.plot$features.plot], levels = unique(feature.groups))
  }
  
  
  data.plot1 <- subset(data.plot,data.plot$id=='CNE')
  data.plot2 <- subset(data.plot,data.plot$id=='Prol1')
  data.plot3 <- subset(data.plot,data.plot$id=='Prol2')
  data.plot4 <- subset(data.plot,data.plot$id=='EMT')
  
  
  print(0)
  data.plot1 <- merge(data.plot1, res.CaCTS1[, c("Name", "LogValue")], 
                      by.x = "features.plot", by.y = "Name", 
                      all.x = TRUE)
  
  
  print(1)
  data.plot2 <- merge(data.plot2, res.CaCTS2[, c("Name", "LogValue")], 
                      by.x = "features.plot", by.y = "Name", 
                      all.x = TRUE)
  
  
  print(2)
  data.plot3 <- merge(data.plot3, res.CaCTS3[, c("Name", "LogValue")], 
                      by.x = "features.plot", by.y = "Name", 
                      all.x = TRUE)
  
  
  print(3)
  data.plot4 <- merge(data.plot4, res.CaCTS4[, c("Name", "LogValue")], 
                      by.x = "features.plot", by.y = "Name", 
                      all.x = TRUE)
  
  
  
  print(4)
  data.plot <- rbind(data.plot1,data.plot2,data.plot3,data.plot4)
  
  print(head(data.plot))
  
  data.plot$pct.exp <- data.plot$LogValue
  # 创建ggplot对象
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) +
    geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    guides(size = guide_legend(title = "Percent Expressed")) +
    labs(x = "Features", y = ifelse(is.null(split.by), "Identity", "Split Identity")) +
    theme_cowplot()
  
  # 添加分面和颜色设置
  if (!is.null(feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(1, "lines"), strip.background = element_blank())
  }
  
  if (split.colors) {
    plot <- plot + scale_color_identity()
  } else if (length(cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  } else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  
  return(plot)
}

 ##Dotplot              
custom_dotplot(seurat.harmony, features = rev(sorted_genes),group.by = 'sampletype',
               assay = 'RNA',scale.min=0.01,dot.limits=c(0.01,75)) + 
  coord_flip() +  # 翻转
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(color = "black", fill = NA, size = 1),  # 添加边框
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +  # 轴标签
  labs(x = NULL, y = NULL) + 
  guides(size = guide_legend("Percent Expression")) +  # legend
  scale_color_gradientn(colours = c("white", "#ce5c69"),limits = c(-2, 2.5))+
  guides(size = guide_legend(title = "Percent Expression")) # 颜色

