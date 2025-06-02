#########################################
##################Figure 8A#############
#########################################
library(tidyverse)
library(ggrepel)

rna <- read.csv('NEC-NE_EMT-NE_limma.csv')
scrna <- read.csv('output_memento.csv')
scrna <- subset(scrna,scrna$subtype=='VR_O')

colnames(scrna)[1] <- 'Gene_symbol'
scrna <- scrna |> 
  mutate(group = case_when( de_coef> 0.2 & de_pval< 0.05  ~ "NEC-dependent (UP)",
                            de_coef < -0.2 & de_pval < 0.05 ~ "EMT-dependent (DN)",
                            TRUE ~ "N.S."))

data <- scrna |> 
  select(Gene_symbol, de_coef, group) |> 
  left_join(rna |> 
              select(Gene_symbol, logFC, adj.P.Val), 
            by = join_by(Gene_symbol)) |> #将RNA和蛋白的FC数据整合到一起
  arrange(de_coef) |> 
  mutate(rank = row_number()) |> #按照RNA数据ERKi_24_hr.logFC排序后得到rank
  filter(adj.P.Val < 0.05) |>
  filter(logFC > 0 | logFC < -0)

data$group <-  factor(data$group, 
                      levels = c("NEC-dependent (UP)",
                                 "EMT-dependent (DN)",
                                 "N.S."))



label <- c("AMIGO2", "ONECUT2","ID1",'DLL3','CHGB','AKR7A3')

label <-intersect(NEP100,data$Gene_symbol)



colors <- c("NEC-dependent (UP)" = "#ce5c69",
            "EMT-dependent (DN)" = "#5770a6",  
            "N.S." = "#a281b1")

ggplot()+
  #绘制所有数据点，x轴是RNA的FC的排序，y轴是蛋白log2FC
  geom_point(data, 
             mapping = aes(x = rank, 
                           y = logFC, 
                           color = group, 
                           fill = group), 
             key_glyph='rect', #指定legend.key的形状为矩形
             shape = 21, 
             size = 1.8,
             stroke = 1.2,
             alpha = 0.8
  ) +
  #绘制APC/C复合体组分基因（带有边缘实线、突出强调）
  geom_point(data |> filter(Gene_symbol %in% NEP100),
             mapping = aes(x = rank, 
                           y = logFC, 
                           fill = group), 
             color = "black",
             shape = 21, 
             size = 1.8,
             stroke = 1.2,
             alpha = 1,
             show.legend = F #不展示legend
  ) +
  #绘制CCNB1，CCNB2和PTTG1的文本标签
  geom_text_repel(data |> filter(Gene_symbol %in% label),
                  mapping  = aes(x = rank, 
                                 y = logFC,
                                 label = Gene_symbol),
                  color = "black", 
                  size = 4.5,
                  min.segment.length = unit(0, "lines"), #定义画出连接线段的最短长短，低于这个长度将不画线段
                  box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.2, "lines"), 
                  segment.color = "black",
                  direction = "both") + #同时在x和y轴方向调整文本标签的方向
  #绘制y = 0的水平方向的虚线
  geom_hline(yintercept = 1, 
             linewidth = 1, 
             linetype = "dotted", 
             color = "grey20") +
  geom_hline(yintercept = -1, 
             linewidth = 1, 
             linetype = "dotted", 
             color = "grey20")+
  #添加APC/C复合体组分基因的图例到图形的右下部分
  annotate("text",                       
           x = 7000,        
           y = -1.5,                         
           label = "NEP100 genes", 
           lineheight = 0.8,
           size = 5,                       
           color = "black",
           hjust = 0) +
  annotate("point",                    
           x = 9600,        
           y = -1.52,   
           shape = 21,
           size = 3,  
           stroke = 1.5,
           fill = "grey60",
           color = "black") +
  #修改坐标轴标题的文本标签
  labs(x = "scRNA NEC_O vs EMT_R (rank of FC)",
       y = "bulk NEC-NE vs EMT-NE (log2FC)", 
       title = "") + 
  scale_y_continuous(limits = c(-2, 2.5), 
                     breaks = c(-2, -1, 1, 2), # 只显示特定的刻度
                     labels = function(x) ifelse(x > -1 & x < 1, "", x)) + # 中间区间不显示标签
  #过滤完后只剩下4000出头的基因，但排序却有10000多个，直接修改标签以和原图一致，我怀疑文章作者也是这么干的！
  #scale_x_continuous(breaks = c(0, 3000, 6000, 9000, 12000),
  #                   labels = c(0, 1000, 2000, 3000, 4000)) +
  scale_color_manual(name = "scRNA",values = colors) +
  scale_fill_manual(name = "scRNA",values = colors) +
  #设置图例中的透明度为1，并覆盖原来图层中的映射（原图层的透明度为0.06）
  guides(fill = guide_legend(override.aes = list(alpha = 1))) + 
  theme_classic() +
  theme(legend.position = c(0.25, 0.95), #将图例置于图形panel中
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size = 14, hjust = 0.5), #图例标题居中
        legend.text = element_text(size = 14),
        legend.key.spacing.y =  unit(0.1, "cm"), #设置图例key在y轴方向的间距
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.5, "cm"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 16)) #设置图标题居中
ggsave("FC-FC dot plot.pdf", height = 4, width = 6,device = cairo_pdf)



#########################################
##################Figure 8C#############
#########################################
library(ggsignif)
df.list <- list()
for (i in 1:5) {
  group <- list_train_vali_Data[[i]][,colnames(list_train_vali_Data[[i]])%in%c('ID','Var')]
  expr <- as.data.frame(t(data.list[[i]]))
  expr$ID <- row.names(expr)
  expr$ID <- gsub("-", ".", expr$ID)
  expr$ID<- gsub("_", ".", expr$ID)
  group$ID <- gsub("-", ".", group$ID)
  group$ID <- gsub("_", ".", group$ID)
  
  expr <- left_join(expr,SU2C)
  
  expr <- expr[,colnames(expr)%in%c('Var','AMIGO2')]
  expr$dataset <- names(list_train_vali_Data)[i]
  
  df.list[[i]] <- expr
  names(df.list)[i] <- names(list_train_vali_Data)[i]
}
names(data.list)
MCTP <- read.csv("G:/importance/undergraduated/数据集整理/bulk data（clinical）/ARPC/CamCap/cli.csv")
expr <- as.data.frame(t(data.list[[9]]))
expr$ID <- row.names(expr)
expr <- left_join(MCTP,expr)
expr <- expr[,colnames(expr)%in%c('group','AMIGO2')]
expr$dataset <- 'CamCap'
expr

list_train_vali_Data <- read_rds("G:/importance/undergraduated/数据集整理/bulk data（clinical）/list_train_vali_Data.rds")

expr <- as.data.frame(t(expr))
expr$ID <- rownames(expr)
group <- list_train_vali_Data[[7]][,c(1,2)]
expr <- left_join(group,expr)

expr <- expr[,colnames(expr)%in%c('Var','AMIGO2')]
names(list_train_vali_Data)[7]
expr$dataset <- 'Pseudobulk'
colnames(expr)[1] <- 'group'

df.combi <- rbind(df.combi,expr)
df.combi <- subset(df.combi,df.combi$dataset%in%c('SU2C-2019','WCDT-MCRPC','WCM','UWRA-CRPC','MDA'))
df.combi$dataset <- factor(df.combi$dataset,levels = c('SU2C-2019','PCaProfilter','WCDT-MCRPC','UWRA-CRPC','WCM','MDA','Pseudobulk'))

write.csv(df.combi,'bulk_AMIGO2.csv')

df.combi <- read.csv('bulk_AMIGO2.csv')
# 创建 KACsig_inPreNeoplasia 图形
ggplot(df.combi, aes(x =group, y = AMIGO2)) +
  geom_boxplot(outlier.size = 1.5, size = 1.2, aes(color = group),outlier.shape = NA,width = 0.6) +  # 设置箱式图的边框颜色
  #geom_line(aes(group = sample), color = '#a281b1', lwd = 1) +  # 添加连接线
  #geom_jitter(aes(fill = group), position = position_jitter(width = 0.3), size = 2, shape = 21, color = "black") +  # 设置散点的边框颜色为黑色
  theme_minimal(base_size = 16, base_family = "sans") +
  theme(
    text = element_text(size = 16),  # 调整字体大小
    axis.title = element_text(size = 24,color = "black"),  # 调整轴标题的大小
    axis.text = element_text(size = 16,color = "black"),  # 调整坐标轴标注文字的大小
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(size = 1),  # 添加坐标轴线
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "",
       y = "AMIGO2 expression(normalization)") +
  scale_color_manual(values =col6, name = NULL)  + # 设置箱式图的边框颜色
  scale_fill_manual(values =col6,name = NULL)  + # 设置散点的填充色
  geom_signif(comparisons = list(c('NEPC','ARPC')), test = "wilcox.test", map_signif_level = F,
              y_position = max(df.combi$AMIGO2))+  # 调整 y_position 以避免与箱式图重叠
  facet_wrap(~ dataset,nrow=1)  # 根据 Group 列分面
col6 <- c('#5770a6','#ce5c69')

#########################################
##################Figure 8D#############
#########################################
library(dplyr)
library(correlation) #install.packages('correlation')
library(ggplot2)
library(ggpubr)
library(viridis)
library(corrplot) #install.packages('corrplot')
library(cols4all)
feature <- c(#AR
  'AR','KLK2','KLK3','TMPRSS2','NKX3-1','NKX3.1','PMEPA1','ALDH1A3','SPDEF',
  #NE 
  'CHGA','CHGB','SYP','ENO2','CHRNB2','ELAVL4','ASCL1','ONECUT2',
  'AMIGO2')
rs <- data.list[[1]][row.names(data.list[[1]])%in%feature,]
rs <- as.data.frame(t(rs))
av.exp <- cor(rs, method= "spearman")
#av.exp
#Heatmap(av.exp,col=colorRampPalette(c("white", "#ce5c69"), alpha = TRUE)(100))

cor <- as.data.frame(av.exp[,'AMIGO2'])
colnames(cor) <- names(data.list)[1]
cor$gene <- row.names(cor)
cor$gene <- as.character(cor$gene)

for (i in 2:18) {
  print(i)
  rs <- data.list[[i]][row.names(data.list[[i]])%in%feature,]
  rs <- as.data.frame(t(rs))
  av.exp <- cor(rs, method= "spearman")
  cor2 <- as.data.frame(av.exp[,'AMIGO2'])
  colnames(cor2) <- names(data.list)[i]
  cor2$gene <- row.names(cor2)
  cor2$gene <- as.character(cor2$gene)
  cor <- left_join(cor,cor2)
}
mycol <- colorRampPalette(c("#5770a6", "white", "#ce5c69"), alpha = TRUE)

row.names(cor) <- cor$gene
cor$gene <- NULL
cor <- as.matrix(cor)

corrplot(cor, method = c('pie'),col=mycol(100),rect.col='black',rect.lwd=3)

row.names(cor) <- cor$gene
cor$gene <- NULL
write.csv(cor,'AMIGO2_cor_final.csv')


#########################################
##################Figure 8F#############
#########################################
###Stacked bar plot###
#Reference:
#Muñoz, K.A., Ulrich, R.J., Vasan, A.K. et al. A Gram-negative-selective antibiotic that spares the gut microbiome. Nature 630, 429–436 (2024).
#https://doi.org/10.1038/s41586-024-07502-0
#Original figure：https://www.nature.com/articles/s41586-024-07502-0/figures/5

#install.packages("pacman") #如果没有pacman，请先安装，pacman是一个管理R包的工具
library(dplyr)
library(tidyverse)
p_load(reshape2)
getwd()
#由于作者没有提供原始数据，这里生成10行，24列的随机数据用于画图
set.seed(123) #设置随机种子，方便复现
dt<- read.csv('IHC.csv')

dt_long <- dt %>%
  pivot_longer(
    cols = c(High.Positive, Positive, Low.Positive, Negative),  # 需要转置的列
    names_to = "IHC",  # 新列，用于存储原列的名称
    values_to = "Value"  # 新列，用于存储原列的值
  )

dt_long <- subset(dt_long,dt_long$Gene=='AR')
dt_long$Value <- dt_long$Value/100
dt_long$IHC <- factor(dt_long$IHC,levels = rev(c('Negative', 'Low.Positive', 'Positive', 'High.Positive')))

col6 <- c('#5770a6','#a281b1','#e0a980','#ce5c69')
ggplot(data = dt_long,aes(x = ID,y = Value, fill = IHC))+
  geom_bar(stat = "identity",position = "stack")+       
  labs(x = "", y = "IHC Relative abundance",fill = "AR expression")+
  scale_fill_manual(values = c("Negative"="#5770a6","Low.Positive"="#a281b1" , "Positive"="#e0a980", "High.Positive"="#ce5c69")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0))+ #limits设置刻度范围；expand = c(lower, upper)：用于设置坐标轴的扩展范围。
  theme_classic()+
  theme(
    plot.subtitle = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 14),
    axis.title =  element_text(color = "black", size = 14),
    axis.text.x = element_blank(), #隐藏X轴文本
    axis.ticks.x = element_blank(),#隐藏X轴刻度
    strip.text.x = element_text(size = 14),
    strip.background = element_blank(),
    strip.placement = "outside",   #facet标签将显示在绘图区域外部
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.4, "cm"),
    legend.key.spacing.y = unit(0.2, "cm"),
    panel.spacing = unit(0.05, "lines"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())+
  facet_grid(~ ID,scales = "free", switch="both")
#分面标签通常展示在plot的上面或右面，switch参数设置为x，上面的标签会放到下面；switch参数设置为y，右面的标签会放到左边，设置为both，相当于x+y
#https://stackoverflow.com/questions/3261597/change-the-position-of-the-strip-label-in-ggplot-from-the-top-to-the-bottom
p
ggsave(filename = "Stacked bar plot.pdf", plot =p, height = 4, width = 6 )

#########################################
##################Figure 8H#############
#########################################
cor_ihc <- read.csv('cor_AR_AMIGO2.csv',row.names = 1)
library(circlize)
col_fun = colorRamp2(c(15,19,23), c('#5770a6','white',"#ce5c69"))
Heatmap(cor_ihc[,-3],  # 不包含 NE 列
        name = "IHC-score",
        col = col_fun,
        na_col = "grey",
        border = TRUE,
        border_gp = gpar(lty = 1, lwd = 2, col = "black"),
        rect_gp = gpar(col = "white"),
        cluster_rows = T,  # 启用行聚类
        cluster_columns = FALSE,  # 禁用列聚类
        row_names_side = "left",
        row_dend_side = "right",
        row_split = cor_ihc$sampletype,
        row_title_side = "left",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "",  # 添加列标题
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_gap = unit(5, "mm"),  # 调整行分割的宽度
        column_gap = unit(3, "mm"),
        heatmap_legend_param = list(title = "IHC-score", title_gp = gpar(fontsize = 12, fontface = "bold")))


#########################################
##################Figure 8I#############
#########################################
#see the code’4.Survival analysis.R‘
#see the code’4.Survival analysis.R‘  
cli <- read.csv('cli.csv')

#########################################
##################Figure 8K#############
#########################################

drug <- read.csv('drug.csv')
drug<- subset(drug,cor>0&fdr<0.05)
CTRP<- subset(drug,dataset=='CTRP')
GDSC<- subset(drug,dataset=='GDSC')

drug_com <- intersect(CTRP$drug,GDSC$drug)

drug2 <- subset(drug,drug%in%drug_com)
drug2$fdr[drug2$fdr == 0] <- 2.887671361E-34

drug$drug <- factor(drug2$drug,levels = c('ZSTK474','PIK-93','PI-103','OSI-027','AZD8055','PHA-793887','AZD7762','UNC0638','TPCA-1','SNX-2112','PAC-1','OSI-930'))
ggplot(drug, aes(x = dataset, y = drug, size = -log2(fdr), color = cor)) +
  geom_point(alpha = 1) +  # 设置气泡透明度
  scale_size_continuous(range = c(2, 8)) +  # 设置气泡大小范围
  scale_color_gradientn(colors = c('white', '#ce5c69'),limits=c(0,0.4)) +  # 设置颜色渐变
  theme_minimal() +  # 使用简约主题
  theme(
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 2),  # 添加边框
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 18,color = "black"),  # 修改 x 轴文本字体大小
    axis.text.y = element_text(size = 18,color = "black"),  # 修改 y 轴文本字体大小
    axis.title.x = element_text(size = 20, face = "bold",color = "black"),  # 修改 x 轴标题字体大小
    axis.title.y = element_text(size = 20, face = "bold",color = "black"),  # 修改 y 轴标题字体大小
    panel.spacing = unit(2, "lines"), 
    strip.background = element_blank()
  ) +
  labs(
    title = "",
    x = "Dataset",
    y = "Drug",
    size = "-log2(FDR)",  # 更新 size 图例标签
    color = "Correlation"  # 更新 color 图例标签
  ) +
  guides(
    size = guide_legend(title = "-log2(FDR)"),  # 修改 size 图例标题
    color = guide_colorbar(title = "Correlation")  # 修改 color 图例标题
  )

write.csv(drug2,'drug2.csv')


# 使用 reshape2 包的 melt 函数将数据从宽格式转换为长格式
data <- read.csv('drug3.csv')

library(reshape2)

melted_data <- melt(data, id.vars = c("drug"))

# 将值转换为数字类型
melted_data$value <- as.numeric(as.character(melted_data$value))

# 提取数据中的唯一变量名
unique_names <- unique(melted_data$drug)
c("ARPC" = "#5770a6",
  "CRPC" = "#a281b1",
  "NEPC" = "#ce5c69")
# 定义更新的颜色向量
colors <- c(
  "#ce5c69", "#5770a6", "#ce5c69", "#5770a6", "#ce5c69", 
  "#5770a6", "#a281b1", "#5770a6", "#5770a6", "#ce5c69", 
  "#ce5c69", "#a281b1"
)

# 创建颜色映射字典
color_dict <- setNames(colors, unique_names)

# 右图：组织2数据的条形图，并在柱子顶端添加圆圈
ggplot(melted_data, 
       aes(x = factor(drug, levels = unique_names), y = value)) +
  geom_bar(stat = "identity", width = 0.1, aes(fill = drug), show.legend = FALSE) +  # 根据 protein 分组着色
  geom_point(aes(y = value, size = value, fill = drug, color = drug), shape = 21) +  # 添加圆圈
  scale_fill_manual(values = color_dict) +  # 设置圆圈填充颜色
  scale_color_manual(values = color_dict) +  # 设置圆圈边框颜色
  scale_size(range = c(2, 8)) +  # 调整圆圈大小范围
  coord_flip() +  # 翻转坐标轴
  theme_minimal() +  # 使用简约主题
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转 x 轴文本
    panel.grid = element_blank(),  # 移除网格线
    axis.text = element_text(color = "black")  # 设置轴文本颜色
  ) +
  labs(
    title = NULL, 
    subtitle = "level_in_tissueB", 
    x = NULL, 
    y = NULL, 
    size = "value of tissueB", 
    color = element_blank(), 
    fill = element_blank()  # 省略坐标轴标签
  )
