###################数据处理#########################
library(dplyr)
library(survival)
library(plyr)
library(dplyr)
library(pROC)
library(survminer)
library(ggsurvfit)
library(survival)
library(ggsurvfit)
library(patchwork)
library(ggpp)

prog.data.list <- readRDS('prog.data.list.rds')
df <- prog.data.list[[1]]

expr<- df[,-c(1:3)]
new <- NEP100(expr,type='bulk',layer='data',species = 'homo')
new$ID <- row.names(new)

cli <- df[,1:3]
cli <- cli[,colnames(cli)%in%c('ID','OS.time','OS')]
colnames(cli)  <- c('ID','OS','OS.time')
cli$ID <- gsub("-", ".", cli$ID)

cli <- subset(cli,cli$ID%in% colnames(new))
row.names(cli) <- cli[,1]
cli <- na.omit(cli)

NEP <- new[,colnames(new)%in% c('ID','NEP100')]

merged_data <- left_join(NEP,cli,by='ID')

df <- merged_data[,colnames(merged_data)%in% c('ID','OS.time','OS','NEP100')]
df <-  na.omit(df)

#df$days_to_recurrence <- df$days_to_recurrence/30
#df <- na.omit(df)
df$OS.time <- as.numeric(df$OS.time)
df$OS <- as.numeric(df$OS)
res.cut <- surv_cutpoint(df,time = "OS.time", event = "OS",variables = 'NEP100')
cutoff <- res.cut$cutpoint$cutpoint
df$group <- ifelse(df$VR.NE >=1, "High", "Low") %>% factor(levels = c("Low","High"))
#########################################
##################Figure 5A##############
#########################################
###########premodel
res <- EZGen::create_model(prog.data.list=prog.data.list,
                           train.data.pos=1,
                           gene.list=NEP100,
                           unicox_pcutoff=0.1,
                           top.num = 100,
                           method='out_RSF',#c('all','out_RSF')
                           hm.col=c("#5770A6", "#FFFFFF", "#CE5C69"),
                           cohort.col=c('#b30c2a','#ce5c69','#e0a980','#f4c889','#bdd5a3','#519981','#8ba1c6','#5770a6','#a281b1'))
write.csv(res,'out_Cindex.csv')

###########presig
folder_path <- './data'
file_list <- list.files(path = folder_path, pattern = "\\.csv$", full.names = T)
file_ID<- list.files(path = folder_path, pattern = "\\.csv$", full.names = F)
c_index.list <- list()
for (i in seq_along(file_list)) {
  df <- read.csv(file_list[i],row.names = 1)
  fit <- coxph(Surv(OS.time,OS)~NEP100,data=df)
  sum.surv <- summary(fit)
  c_index <-sum.surv$concordance
  c_index.list[[i]] <- c_index
  names(c_index.list)[i] <- file_ID[i]
}
saveRDS(c_index.list,'presig.rds')
#########################################
##################Figure 5B##############
#########################################
fit <- survfit(Surv(OS.time, OS)~group,
               data= df)
ggsurvplot(fit, conf.int=F, pval=TRUE,risk.table = F,)


unicox <- function(df=df,ncol=5){
  cont_covariates <- colnames(df)[c()]
  cate_covariates <- colnames(df)[c(ncol)]
  
  Uni_cox_model<- function(x){
    surv <- as.formula(paste0 ("y~",x))
    cox <- coxph(surv,data=df)
    sum <- summary(cox)
    HR <- round(sum$coefficients[,2],3)
    PValue <- sum$coefficients[, 5]  # 保留原始P值
    PValue_formatted <- format(PValue, nsmall = 3)  # 确保显示三位小数
    lower <- round(sum$conf.int[,3],3)
    upper <- round(sum$conf.int[,4],3)
    subchar <- rownames(sum$coefficients)
    HRa <- paste0(sprintf("%1.3f",HR) ," [",sprintf("%1.3f",lower),", ",sprintf("%1.3f",upper),"]")
    Uni_cox_model<- data.frame('Characteristics' = paste(x,"_",subchar,sep=""),
                               'HRa'=HRa,
                               'PValue' = PValue,
                               'HR' = HR,
                               'lower' = lower,
                               'upper' = upper
    )
    return(Uni_cox_model)
  }
  
  #转换成数据框，并转置
  univ_results  <- do.call(rbind,lapply(c(cont_covariates,cate_covariates) %>% na.omit ,Uni_cox_model))
  #最后，将P值=0的变为p<0.0001
  univ_results$PValue <- as.numeric(as.character(univ_results$PValue))
  names(univ_results) <- c("Variants","Hazard Ratio (95%CI) ","P-value ","HR ","Lower ","Upper ")
  return(univ_results)
}
y<- Surv(time=df$OS.time,event=df$OS==1)#1为事件发生
Camcap.COX <- unicox(df=df,ncol=6)
Camcap.COX

unicox.rs.res <- read.csv('ALL_cox.csv',row.names = 1)
metamodel <- Mime1::cal_unicox_meta_ml_res(input = unicox.rs.res)
Mime1::meta_unicox_vis(metamodel,
                       dataset = names(prog.data.list))
#########################################
##################Figure 5C##############
#########################################

max(df$OS.time)
survfit2(Surv(OS.time, OS) ~ group, data = df) %>%
  ggsurvfit(linewidth = 2) +
  add_pvalue(
    location = "annotation",
    caption = "{p.value}",
    prepend_p = TRUE,
    pvalue_fun = format_p,
    rho = 0
  )+
  add_risktable(
    risktable_height = 0.15,
    risktable_stats = c("{n.risk}"),
    stats_label = list(n.risk = "No. at risk"),
    size = 5,
    theme = list(
      theme_risktable_default(axis.text.y.size = 14, plot.title.size = 14),
      theme(plot.title = element_text(face = "plain"))
    )) +
  add_risktable_strata_symbol(symbol = "\U25CF", size = 20) +
  add_censor_mark(size = 5, shape = 73) +
  #add_quantile(y_value = 1, linetype = "dashed", color = "black", linewidth = 2) +
  labs(title = "",
       x = "Time (months)",
       y = "OS (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 60, 15), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25),expand = c(0, 0.01)) +
  scale_color_manual(values = c('#5770a6', '#ce5c69')) +
  scale_fill_manual(values = c('#5770a6', '#ce5c69')) +
  guides(color = guide_legend(ncol = 1)) +
  theme_classic() +
  theme(axis.text = element_text(size = 20,color = "black"),
        axis.ticks.length = unit(2, "mm"),
        axis.title.y = element_text(size = 24,color = "black"),
        axis.title.x = element_text(size = 24,color = "black"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 13,color = "black"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", color = "grey"),  # 添加主网格线
        panel.grid.minor = element_line(size = 0.25, linetype = "dotted", color = "lightgrey"),  # 添加次网格线
        legend.background = element_blank(),
        legend.position = c(0.8,0.8),
        panel.border = element_rect(size = 2, color = "black", fill = NA))
ggsave(file = "Nature Medicine.pdf", width = 10, height = 6, device = cairo_pdf)

#########################################
##################Figure D-E##############
#########################################
#see the code’1.Single-cell atlas.R‘
#see the code’1.Single-cell atlas.R‘    
#########################################
##################Figure 5F##############
#########################################
library(IOBR)
library(correlation) #install.packages('correlation')
library(ggplot2)
library(ggpubr)
library(ggpointdensity) #install.packages('ggpointdensity')
library(viridis)
library(corrplot) #install.packages('corrplot')
library(corrmorant) #remotes::install_github('r-link/corrmorant')
library(cols4all)
library(ComplexHeatmap)
library(circlize)

nep <- NEP100(data.list[[1]])
rs <- IOBR::calculate_sig_score(pdata           = NULL,
                                eset            = data.list[[1]],
                                signature       = c(pathway_list),
                                method          = "ssGSEA",
                                mini_gene_count = 0,
                                adjust_eset = T,
                                parallel.size = 8)

nep$ID <- rownames(nep)
nep <- nep[,c('NEP100','ID')]
rs <- left_join(rs,nep,by='ID')
rs <- as.data.frame(rs)
row.names(rs) <- rs[,1]
rs <- rs[,-1]
av.exp <- cor(rs, method= "spearman")

cor <- as.data.frame(av.exp[,'NEP100'])
colnames(cor) <- names(data.list)[1]
cor$pathway <- row.names(cor)

for (i in 2:18) {
  print(i)
  nep <- NEP100(data.list[[i]])
  rs <- IOBR::calculate_sig_score(pdata           = NULL,
                                  eset            = data.list[[i]],
                                  signature       = c(pathway_list),
                                  method          = "ssGSEA",
                                  mini_gene_count = 0,
                                  adjust_eset = T,
                                  parallel.size = 8)
  
  nep$ID <- rownames(nep)
  nep <- nep[,c('NEP100','ID')]
  rs <- left_join(rs,nep,by='ID')
  rs <- as.data.frame(rs)
  row.names(rs) <- rs[,1]
  rs <- rs[,-1]
  av.exp <- cor(rs, method= "spearman")
  cor2 <- as.data.frame(av.exp[,'NEP100'])
  colnames(cor2) <- names(data.list)[i]
  cor2$pathway <- row.names(cor2)
  cor2$pathway <- as.character(cor2$pathway)
  cor <- left_join(cor,cor2,by='pathway')
}
rownames(cor) <- cor$pathway
cor$pathway <- NULL
write.csv(cor,'cor.csv')

Heatmap(cor,  # 不包含 NE 列
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
        #row_split = row_split,
        #column_split = column_clusters2$x,
        #row_title_side = "right",
        row_title_gp = gpar(fontsize = 12, fontface = "bold"),
        show_row_names = T,
        show_column_names = T,
        column_title = "",  # 添加列标题
        column_title_side = "top",
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        row_gap = unit(10, "mm"),  # 调整行分割的宽度
        column_gap = unit(3, "mm"))  # 调整列分割的宽度
# right_annotation = ha_right,
#top_annotation = ha_top_combined)

# 创建pie热图
mycol <- colorRampPalette(c("#5770a6", "white", "#ce5c69"), alpha = TRUE)
corrplot(cor, method = c('pie'),col=mycol(100))

