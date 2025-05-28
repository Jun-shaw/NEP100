##################建模####################
##############函数
cal_sig <- function (train_data, list_train_vali_Data, candidate_genes = NULL, 
                     methods = NULL, seed = 5201314, cores_for_parallel = 12) 
{
  message("---loading the packages ---")
  if (T) {
    library(stringr)
    library(gridExtra)
    library(future)
    library(sva)
    library(e1071)
    library(pROC)
    library(ROCit)
    library(caret)
    library(doParallel)
    library(cancerclass)
    library(dplyr)
  }
  model.Dev <- function(training, method, sig) {
    training <- training[, colnames(training) %in% c("Var", 
                                                     sig)]
    Grid <- list(nb = expand.grid(fL = c(0, 0.5, 1, 1.5, 
                                         2), usekernel = TRUE, adjust = c(0.5, 0.75, 1, 1.25, 
                                                                          1.5)), svmRadialWeights = expand.grid(sigma = c(5e-04, 
                                                                                                                          0.001, 0.005, 0.01, 0.05), C = c(1, 3, 5, 10, 20), 
                                                                                                                Weight = c(0.1, 0.5, 1, 2, 3, 5, 10)), rf = expand.grid(mtry = c(2, 
                                                                                                                                                                                 42, 83, 124, 165, 205, 246, 287, 328, 369)), kknn = expand.grid(kmax = c(5, 
                                                                                                                                                                                                                                                          7, 9, 11, 13), distance = 2, kernel = "optimal"), 
                 adaboost = expand.grid(nIter = c(50, 100, 150, 200, 
                                                  250), method = c("Adaboost.M1", "Real adaboost")), 
                 LogitBoost = expand.grid(nIter = c(11, 21, 31, 41, 
                                                    51, 61, 71, 81, 91, 101)))
    TuneLength <- list(nb = nrow(Grid[["nb"]]), svmRadialWeights = nrow(Grid[["svmRadialWeights"]]), 
                       rf = nrow(Grid[["rf"]]), kknn = nrow(Grid[["kknn"]]), 
                       adaboost = nrow(Grid[["adaboost"]]), LogitBoost = nrow(Grid[["LogitBoost"]]))
    ls_model <- lapply(method, function(m) {
      if (m == "cancerclass") {
        pData <- data.frame(class = training$Var, sample = rownames(training), 
                            row.names = rownames(training))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(training[, -1])
        Sig.Exp.train <- ExpressionSet(assayData = as.matrix(Sig.Exp), 
                                       phenoData = phenoData)
        predictor <- fit(Sig.Exp.train, method = "welch.test")
        model.tune <- predictor
      }
      else {
        f <- 5
        r <- 10
        n <- f * r
        seeds <- vector(mode = "list", length = n + 1)
        for (i in 1:n) seeds[[i]] <- sample.int(n = 1000, 
                                                TuneLength[[m]])
        seeds[[n + 1]] <- sample.int(1000, 1)
        ctrl <- trainControl(method = "repeatedcv", number = f, 
                             summaryFunction = twoClassSummary, classProbs = TRUE, 
                             repeats = r, seeds = seeds)
        model.tune <- train(Var ~ ., data = training, 
                            method = m, metric = "ROC", trControl = ctrl, 
                            tuneGrid = Grid[[m]])
      }
      print(m)
      return(model.tune)
    })
    names(ls_model) <- method
    return(ls_model)
  }
  cal.model.auc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% 
                                   c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", 
                                                        "Y"))
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    auc <- lapply(1:length(models), function(i) {
      if (models[i] == "cancerclass") {
        model.tune <- ls_model[[i]]
        pData <- data.frame(class = validation$Var, sample = rownames(validation), 
                            row.names = rownames(validation))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- ExpressionSet(assayData = as.matrix(Sig.Exp), 
                                      phenoData = phenoData)
        prediction <- predict(model.tune, Sig.Exp.test, 
                              "N", ngenes = nrow(Sig.Exp), dist = "cor")
        roc <- roc(response = prediction@prediction[, 
                                                    "class_membership"], predictor = as.numeric(prediction@prediction[, 
                                                                                                                      "z"]))
        roc_result <- coords(roc, "best")
        auc <- data.frame(ROC = as.numeric(roc$auc), 
                          Sens = roc_result$sensitivity[1], Spec = roc_result$specificity[1])
      }
      else {
        model.tune <- ls_model[[i]]
        prob <- predict(model.tune, validation[, -1], 
                        type = "prob")
        pre <- predict(model.tune, validation[, -1])
        test_set <- data.frame(obs = validation$Var, 
                               N = prob[, "N"], Y = prob[, "Y"], pred = pre)
        auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
      }
      return(auc)
    }) %>% base::do.call(rbind, .)
    rownames(auc) <- names(ls_model)
    return(auc)
  }
  cal.model.roc <- function(res.by.model.Dev, cohort.for.cal, sig) {
    library(dplyr)
    rownames(cohort.for.cal) <- cohort.for.cal$ID
    validation <- cohort.for.cal[, colnames(cohort.for.cal) %in% 
                                   c("Var", sig)]
    validation$Var <- factor(validation$Var, levels = c("N", 
                                                        "Y"))
    ls_model <- res.by.model.Dev
    models <- names(ls_model)
    roc <- lapply(1:length(models), function(i) {
      if (!models[i] == "cancerclass") {
        prob <- predict(ls_model[[models[i]]], validation[, 
                                                          -1], type = "prob")
        pre <- predict(ls_model[[models[i]]], validation[, 
                                                         -1])
        test_set <- data.frame(obs = validation$Var, 
                               N = prob[, "N"], Y = prob[, "Y"], pred = pre)
        roc <- ROCit::rocit(score = test_set$N, class = test_set$obs, 
                            negref = "Y")
      }
      else {
        pData <- data.frame(class = validation$Var, sample = rownames(validation), 
                            row.names = rownames(validation))
        phenoData <- new("AnnotatedDataFrame", data = pData)
        Sig.Exp <- t(validation[, -1])
        Sig.Exp.test <- ExpressionSet(assayData = as.matrix(Sig.Exp), 
                                      phenoData = phenoData)
        prediction <- predict(ls_model[[models[i]]], 
                              Sig.Exp.test, "N", ngenes = nrow(Sig.Exp), 
                              dist = "cor")
        roc <- roc(response = prediction@prediction[, 
                                                    "class_membership"], predictor = as.numeric(prediction@prediction[, 
                                                                                                                      "z"]))
      }
    })
    names(roc) <- models
    return(roc)
  }
  message("---loading the function---")
  common_feature <- c("ID", "Var", candidate_genes)
  colnames(train_data) <- gsub("-", ".", colnames(train_data))
  candidate_genes <- gsub("-", ".", candidate_genes)
  for (i in names(list_train_vali_Data)) {
    common_feature <- intersect(common_feature, colnames(list_train_vali_Data[[i]]))
  }
  if (all(is.element(methods, c("nb", "svmRadialWeights", "rf", 
                                "kknn", "adaboost", "LogitBoost", "cancerclass")))) {
    list_train_vali_Data <- lapply(list_train_vali_Data, 
                                   function(x) {
                                     x[, common_feature]
                                   })
    list_train_vali_Data <- lapply(list_train_vali_Data, 
                                   function(x) {
                                     x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, as.numeric)
                                     rownames(x) <- x$ID
                                     return(x)
                                   })
    list_train_vali_Data <- lapply(list_train_vali_Data, 
                                   function(x) {
                                     x[, c(1:2)] <- apply(x[, c(1:2)], 2, as.factor)
                                     return(x)
                                   })
    list_train_vali_Data <- lapply(list_train_vali_Data, 
                                   function(x) {
                                     x <- x[!is.na(x$Var) & !is.na(x$Var), ]
                                     return(x)
                                   })
    list_train_vali_Data <- lapply(list_train_vali_Data, 
                                   function(x) {
                                     x[, -c(1:2)] <- apply(x[, -c(1:2)], 2, function(x) {
                                       x[is.na(x)] <- mean(x, na.rm = T)
                                       return(x)
                                     })
                                     return(x)
                                   })
    train_data <- train_data[, common_feature]
    train_data[, -c(1:2)] <- apply(train_data[, -c(1:2)], 
                                   2, as.numeric)
    train_data[, c(1:2)] <- apply(train_data[, c(1:2)], 2, 
                                  as.factor)
    rownames(train_data) <- train_data$ID
    est_dd <- as.data.frame(train_data)[, common_feature[-1]]
    pre_var <- common_feature[-c(1:2)]
    print(paste0("There existing ", length(candidate_genes), 
                 " genes in candidate genes"))
    print("Intersetion of the candidate genes and the colnames of the provided data")
    print(paste0("There existing ", length(pre_var), " genes in candidate genes, colnames of training data, colnames of validation data"))
    cl <- makePSOCKcluster(cores_for_parallel)
    registerDoParallel(cl)
    res.model <- model.Dev(training = train_data, method = methods, sig = pre_var)
    stopCluster(cl)
    ml.auc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.auc(res.by.model.Dev = res.model, 
                               cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.auc) <- names(list_train_vali_Data)
    ml.roc <- lapply(list_train_vali_Data, function(x) {
      res.tmp <- cal.model.roc(res.by.model.Dev = res.model, 
                               cohort.for.cal = x, sig = pre_var)
      return(res.tmp)
    })
    names(ml.roc) <- names(list_train_vali_Data)
    res <- list()
    res[["model"]] <- res.model
    res[["auc"]] <- ml.auc
    res[["roc"]] <- ml.roc
    res[["sig.gene"]] <- pre_var
    return(res)
  }
  else {
    print("Please provide the correct parameters")
  }
}

res.ici <- cal_sig(train_data = list_train_vali_Data[[1]],
                   list_train_vali_Data = list_train_vali_Data,
                   candidate_genes = genelist,
                   methods = c('svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                   seed = 20000709,
                   cores_for_parallel = 4)
auc_vis_category_all(model.res,dataset = names(list_train_vali_Data),
                     order= names(list_train_vali_Data))
saveRDS(res.ici,'model.res.rds')



################ 重要性plot ################
imp <- model.res[["model"]][["rf"]][["finalModel"]][["importance"]]
imp <- as.data.frame(imp)

# 确保列名正确
imp$feature <- row.names(imp)
colnames(imp)[1] <- 'importance'

library(ggplot2)
library(dplyr)

# 假设 intersection_down_sample 和 intersection_up_sample 已经定义
imp <- imp %>%
  mutate(group = case_when(
    feature %in% intersection_down_sample ~ 'down',
    feature %in% intersection_up_sample ~ 'up',
    TRUE ~ NA_character_  # 处理未分类的特征
  ))

# 确保为数值型
imp$importance <- as.numeric(as.character(imp$importance))  # 转换为数值型
imp$feature <- as.factor(imp$feature)                        # 转换为因子型
imp$group <- as.factor(imp$group)                            # 转换为因子型

# 筛选前 20 个 'up' 特征
top_up_features <- imp %>%
  filter(group == 'up') %>%
  arrange(desc(importance)) %>%
  slice_head(n = 20) %>%
  pull(feature)

# 只保留 'up' 和 'down' 特征
imp <- imp %>%
  filter(group == 'up' & feature %in% top_up_features | group == 'down')

# 将 'down' 特征的 importance 值取负值，以便在图的左侧显示
imp <- imp %>%
  mutate(importance = ifelse(group == 'down', -importance, importance))

# 将 importance 转换为百分比并保留两位小数
imp$importance_pct <- round(imp$importance * 100, 2)

ggplot(data = imp, aes(x = importance, y = reorder(feature, importance))) +
  geom_col(aes(fill = importance), size = 0.25, color = "white") +  # 使用 importance 作为填充
  
  # 在条形上显示特征名称和百分比
  geom_text(aes(x = ifelse(group == 'down', 0.02, -0.02), 
                y = feature, 
                label = feature,  # 显示特征名称
                hjust = ifelse(group == 'down', 0, 1)),  # 'down' 组左对齐，其他组右对齐
            vjust = 0.5, size = 4) +
  
  geom_text(aes(x = ifelse(group == 'down', importance * .3, importance * .3), 
                y = feature, 
                label = paste0(importance_pct, "%"),  # 显示特征名称和百分比
                hjust = ifelse(group == 'down', 1, 0)), 
            vjust = .4, size = 4, col = 'white', fontface = "bold") +
  
  geom_text(aes(x = ifelse(group == 'down', 0.01, -0.01), 
                y = feature, 
                label = '--',  # 显示特征名称和百分比
                hjust = ifelse(group == 'down', 1, 0)), 
            vjust = .4, size = 4) +
  
  geom_vline(xintercept = 0, size = 1, color = "black") +
  geom_hline(yintercept = 0, size = 1, color = "black") +
  geom_segment(data = data.frame(x = seq(-0.5, max(imp$importance), by = 0.05), 
                                 y = 0, 
                                 yend = 0.2), size = 1, 
               aes(x = x, xend = x, y = y, yend = yend), 
               color = "black") +  # 短竖线的颜色和高度
  
  scale_y_discrete(expand = c(.025, .025)) +
  
  scale_fill_gradient2(low = "#364888", mid = "white", high = "#b30c2a", midpoint = 0, guide = "none") + # 渐变色
  
  coord_cartesian(clip = "off") +
  
  theme_minimal() +
  
  theme(panel.grid = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(face = "bold", size = rel(1), color = "black"))
#################NEP100###############
NEP100 <- function(data=SU2C,type='bulk',assay="RNA",layer='data',species='homo'){
  library(Seurat)
  library(ggpubr)
  if (species=='homo') {
    NEP100.UP <- c("ACYP1", "ADCYAP1", "AKR7A3", "AMIGO2", "APLP1", "ASF1B", "ASPHD1",
                   "ASPM", "BARD1", "BEX1", "CA8", "CAMK2N2", "CCNE2", "CDH2", "CDK1",
                   "CDKN2A", "CDKN2C", "CDKN3", "CDT1", "CELF4", "CENPF", "CENPK",
                   "CENPW", "CHEK1", "CHGA", "CHGB", "CLSPN", "CRMP1", "CRYBA2",
                   "DLL3", "DNER", "E2F1", "ECT2", "ELAVL3", "ESPL1", "FAM111B",
                   "FBXO5", "FIGNL1", "FOXA2", "GGH", "GINS2", "GPSM2", "HELLS",
                   "HEPACAM2", "HMGB3", "INSM1", "LCN15", "MCM4", "MEST", "MLLT11",
                   "MYBL2", "NOL4", "NPTX1", "NUF2", "NUSAP1", "ORC6", "PEG10",
                   "PRC1", "PROX1", "PSMC3IP", "PTTG1", "RFC4", "RIMKLA", "RIMS2",
                   "RUNDC3A", "SCG2", "SCG3", "SCGN", "SMC2", "SMC4", "SNAP25",
                   "SOX2", "SRRM3", "ST18", "STMN1", "STX1A", "SYP", "TAGLN3",
                   "TMEM176A", "TOP2A", "TPX2", "TRIM9", "TUBB2B", "TYMS", "UBE2C",
                   "UBE2T", "UCHL1", "UHRF1", "VGF", "ZWINT")
    NEP100.DN <- c("KLK3", "KLK4", "MSMB", "RAB27B", "SLC44A4", "SLC45A3", "SPDEF", "TMPRSS2", "TSPAN1", "ZG16B")
  }
  if (species=='mice') {
    NEP100.UP <- c("Adcyap1", "Amigo2", "Aplp1", "Asphd1", "Aspm",
                   "Bard1", "Bex1", "Car8", "Camk2n2", "Ccne2", "Cdh2", "Cdk1",
                   "Cdkn2c", "Cdkn3", "Cdt1", "Celf4", "Cenpf", "Cenpk", "Cenpw", "Chek1",
                   "Chga", "Chgb", "Clspn", "Crmp1", "Cryba2", "Dll3", "Dner", "E2f1",
                   "Ect2", "Elavl3", "Espl1", "Fbxo5", "Fignl1", "Foxa2", "Ggh",
                   "Gins2", "Gpsm2", "Hells", "Hepacam2", "Hmgb3", "Insm1",  "Mcm4",
                   "Mest", "Mllt11", "Mybl2", "Nol4", "Nptx1", "Nuf2", "Nusap1", "Orc6",
                   "Peg10", "Prc1", "Prox1", "Psmc3ip", "Pttg1", "Rfc4", "Rimkla", "Rims2",
                   "Rundc3a", "Scg2", "Scg3", "Scgn", "Smc2", "Smc4", "Snap25", "Sox2",
                   "Srrm3", "St18", "Stmn1", "Stx1a", "Syp", "Tagln3", "Tmem176a", "Top2a",
                   "Tpx2", "Trim9", "Tubb2b", "Tyms", "Ube2c", "Ube2t", "Uchl1", "Uhrf1",
                   "Vgf", "Zwint")
    NEP100.DN <- c("Klk4", "Msmb", "Rab27b", "Slc45a3", "Spdef", "Tmprss2", "Tspan1")
  }
  
  NEP100<- c(NEP100.UP,NEP100.DN)
  imp <- data.frame(
    feature = c("ACYP1", "ADCYAP1", "AKR7A3", "AMIGO2", "APLP1", "ASF1B",
                "ASPHD1", "ASPM", "BARD1", "BEX1", "CA8", "CAMK2N2",
                "CCNE2", "CDH2", "CDK1", "CDKN2A", "CDKN2C", "CDKN3",
                "CDT1", "CELF4", "CENPF", "CENPK", "CENPW", "CHEK1",
                "CHGA", "CHGB", "CLSPN", "CRMP1", "CRYBA2", "DLL3",
                "DNER", "E2F1", "ECT2", "ELAVL3", "ESPL1", "FAM111B",
                "FBXO5", "FIGNL1", "FOXA2", "GGH", "GINS2", "GPSM2",
                "HELLS", "HEPACAM2", "HMGB3", "INSM1", "LCN15", "MCM4",
                "MEST", "MLLT11", "MYBL2", "NOL4", "NPTX1", "NUF2",
                "NUSAP1", "ORC6", "PEG10", "PRC1", "PROX1", "PSMC3IP",
                "PTTG1", "RFC4", "RIMKLA", "RIMS2", "RUNDC3A", "SCG2",
                "SCG3", "SCGN", "SMC2", "SMC4", "SNAP25", "SOX2",
                "SRRM3", "ST18", "STMN1", "STX1A", "SYP", "TAGLN3",
                "TMEM176A", "TOP2A", "TPX2", "TRIM9", "TUBB2B", "TYMS",
                "UBE2C", "UBE2T", "UCHL1", "UHRF1", "VGF", "ZWINT",
                "KLK3", "KLK4", "MSMB", "RAB27B", "SLC44A4", "SLC45A3",
                "SPDEF", "TMPRSS2", "TSPAN1", "ZG16B"),
    importance = c(0.32078655, 0.31467747, 0.22240419, 0.33860124, 0.10814995,
                   0.17895444, 0.11798070, 0.20234881, 0.08742063, 0.11181459,
                   0.42549497, 0.14971475, 0.20739310, 0.36079247, 0.12556780,
                   0.09434643, 0.17682817, 0.23745165, 0.18944805, 0.17280171,
                   0.17339193, 0.16794163, 0.12894988, 0.35231675, 0.22730971,
                   0.12730479, 0.29733918, 0.19465157, 0.13449769, 0.25735914,
                   0.13200176, 0.16006157, 0.11941504, 0.32603263, 0.20330564,
                   0.17097028, 0.20369654, 0.18155342, 0.60661272, 0.19703199,
                   0.13161035, 0.22191526, 0.34497154, 0.20023709, 0.14676781,
                   0.29982674, 0.14644647, 0.41606185, 0.30042446, 0.12290139,
                   0.19216430, 0.19096064, 0.18434198, 0.19350472, 0.15532995,
                   0.12643310, 0.36583652, 0.18277790, 0.44248601, 0.20551648,
                   0.16000738, 0.14208458, 0.22619818, 0.18608543, 0.22547833,
                   0.28218170, 0.15475816, 0.16277019, 0.11682926, 0.15414170,
                   0.16236257, 0.25661811, 0.12413896, 0.24917159, 0.15131578,
                   0.17098293, 0.16769959, 0.19497671, 0.26660330, 0.11860284,
                   0.09370752, 0.35123755, 0.26024785, 0.12681748, 0.13344663,
                   0.24058303, 0.12587054, 0.10406710, 0.12709524, 0.19890655,
                   0.24446303, 0.23863716, 0.16118374, 0.42051391, 0.31661908,
                   0.25063992, 0.34211740, 0.15609867, 0.28086769, 0.19444105)
  )
  
  new <- data[row.names(data) %in%  NEP100,]
  new.up <- new[row.names(new) %in%  NEP100.UP,]
  new.down <- new[row.names(new) %in%  NEP100.DN,]
  if (type=='single') {
    up_data <- t(as.data.frame(GetAssayData(new.up,assay=assay,layer=layer) ) ) # 转换为数值型数据框
    down_data <- t(as.data.frame(GetAssayData(new.down,assay=assay,layer=layer) ) ) # 转换为数值型数据框
    
    imp.up <- imp[1:90,]
    
    new$NEP100.UP <- apply(up_data, 1, function(gene_expression) {
      NEP100.UP <- sum(gene_expression * imp.up$importance, na.rm = TRUE)
      return(NEP100.UP)
    })
    
    imp.down <- imp[91:100,]
    
    new$NEP100.DN <- apply(down_data, 1, function(gene_expression) {
      NEP100.DN <- sum(gene_expression * imp.down$importance, na.rm = TRUE)
      return(NEP100.DN)
    })
    new$NEP100 <- new$NEP100.UP -new$NEP100.DN
    return(new@meta.data)
    
  }
  else if(type=='bulk'){new <- as.data.frame(t(new))
  
  up_data <- as.data.frame(t(new.up))
  down_data <- as.data.frame(t(new.down))
  
  
  
  imp.up <- imp[1:90,]
  
  new$NEP100.UP <- apply(up_data, 1, function(gene_expression) {
    gene_expression <- as.numeric(gene_expression)
    importance_values <- imp.up$importance[match(colnames(up_data), imp.up$feature)]
    NEP100.UP <- sum(gene_expression * importance_values, na.rm = TRUE)
    return(NEP100.UP)
  })
  
  imp.down <- imp[91:100,]
  
  new$NEP100.DN <- apply(down_data, 1, function(gene_expression) {
    gene_expression <- as.numeric(gene_expression)
    importance_values <- imp.down$importance[match(colnames(down_data), imp.down$feature)]
    NEP100.DN <- sum(gene_expression * importance_values, na.rm = TRUE)
    return(NEP100.DN)
  })
  new$NEP100 <- new$NEP100.UP -new$NEP100.DN
  
  return(new)
  
  }
  else{print('type choose single or bulk')}
  
}
###############NEP100跟先前对比#############
library(Mime1)
imp <- read.csv('importance.csv')
list_train_vali_Data <- readRDS('G:/importance/undergraduated/数据集整理/bulk data（clinical）/list_train_vali_Data.rds')
other.res <- cal_other_sig(list_train_vali_Data = list_train_vali_Data,imp=imp)
saveRDS(other.res,'other.res.rds')

data_list <- lapply(other.res$auc.list.up, function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(x[, 1])  # 提取第一列的所有值
  } else {
    stop("The element is not a single-column matrix.")
  }
})
data_list2 <- lapply(other.res$auc.list.all, function(x) {
  if (is.matrix(x) && ncol(x) == 1) {
    return(x[, 1])  # 提取第一列的所有值
  } else {
    stop("The element is not a single-column matrix.")
  }
})

# 合并为矩阵
data_matrix <- do.call(cbind, data_list)
data_matrix2 <- do.call(cbind, data_list2)
# 设置列名

colnames(data_matrix) <- names(other.res$auc.list.up)
colnames(data_matrix2)  <- names(other.res$auc.list.all)
# 查看结果
other.res$auc.this[,1]
data <- cbind(data_matrix,data_matrix2,other.res$auc.this[,1])

print(data)
colnames(data)[17] <- 'This.study'
cindex <- t(data)
cindex[cindex < 0.5] <- 1 - cindex[cindex < 0.5]
cindex <- as.data.frame(cindex)
write.csv(cindex,'other.cindex.csv')

################ AUC热图 ################
cindex <- read.csv('model.cindex.csv',row.names = 1)
cohort.levels <- colnames(cindex)

Cindex_mat=cindex
# 计算每种算法在所有队列中平均C-index
avg_Cindex <- apply(Cindex_mat, 1, mean)
# 对各算法C-index由高到低排序
avg_Cindex <- sort(avg_Cindex, decreasing = T)
# 对C-index矩阵排序
Cindex_mat <- Cindex_mat[names(avg_Cindex), ]
# 保留三位小数
avg_Cindex <- as.numeric(format(avg_Cindex, digits = 3, nsmall = 3))
row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8, border = FALSE,
                                          gp = gpar(fill = "steelblue", col = NA),
                                          add_numbers = T, numbers_offset = unit(-10, "mm"),
                                          axis_param = list("labels_rot" = 0),
                                          numbers_gp = gpar(fontsize = 9, col = "white"),
                                          width = unit(3, "cm")),
                       show_annotation_name = F)


cohort.col <- c('#b30c2a','#ce5c69','#e0a980','#f4c889','#bdd5a3','#5770a6','#a281b1')
cohort.col <- cohort.col[1:length(cohort.levels)]
names(cohort.col) <- colnames(Cindex_mat)
col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                          col = list("Cohort" = cohort.col),
                          show_annotation_name = F)
col_ha@anno_list[["Cohort"]]@color_mapping@levels <- cohort.levels

cellwidth = 1
cellheight = 0.5
library(circlize)
hm.col2 <- colorRamp2(hm.limit, 
                      hm.col)  # 使用蓝色到红色的渐变
hm <- Heatmap(as.matrix(Cindex_mat), name = "C-index",
              right_annotation = row_ha,
              top_annotation = col_ha,
              col = hm.col2,
              rect_gp = gpar(col = "black", lwd = 1), # 边框设置为黑色
              cluster_columns = FALSE, cluster_rows = FALSE, # 不进行聚类，无意义
              show_column_names = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
              height = unit(cellheight * nrow(Cindex_mat), "cm"),
              column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
              column_title = NULL,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path(file.name), width = cellwidth * ncol(Cindex_mat) + 5, height = cellheight * nrow(Cindex_mat) * 0.45)
draw(hm)
invisible(dev.off())



##################验证#############
##########单细胞

new <- NEP100(GSE206962,type='single',layer='data',species = 'homo')
new.PN <- subset(new, subset = celltype %in% c('Luminal', 'NE'))
ggplot(new.PN@meta.data, aes(x =celltype, y = NEP100)) +
  geom_boxplot(outlier.size = 1.5, size = 1.2, aes(color =celltype),outlier.shape = NA,width = 0.6) +  # 设置箱式图的边框颜色
  theme_minimal(base_size = 16, base_family = "sans") + 
  theme(
    text = element_text(size = 16),  # 调整字体大小
    axis.title = element_text(size = 24,color = "black"),  # 调整轴标题的大小
    axis.text = element_text(size = 20,color = "black"),  # 调整坐标轴标注文字的大小
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(size = 1),  # 添加坐标轴线
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "",
       y = "NEP100") +
  scale_color_manual(values =c('#ce5c69','#5770a6'), name = NULL)  + # 设置箱式图的边框颜色
  scale_fill_manual(values = c('#ce5c69','#5770a6'),name = NULL)  + # 设置散点的填充色
  #scale_y_continuous(limits = c(-0.2, 0.3)) +
  geom_signif(comparisons = list(c('NE', 'Luminal')), test = "wilcox.test", map_signif_level = F,
              y_position = max(new.PN@meta.data$NEP100))
#########bulk
new <- NEP100(expr,type='bulk',species = 'homo')
new$ID <- row.names(new)
new.PN <- left_join(new,group,by='ID')

ggplot(new.PN, aes(x =group, y = NEP100)) +
  geom_boxplot(outlier.size = 1.5, size = 1.2, aes(color =group),outlier.shape = NA,width = 0.6) +  # 设置箱式图的边框颜色
  theme_minimal(base_size = 16, base_family = "sans") + 
  theme(
    text = element_text(size = 16),  # 调整字体大小
    axis.title = element_text(size = 24,color = "black"),  # 调整轴标题的大小
    axis.text = element_text(size = 20,color = "black"),  # 调整坐标轴标注文字的大小
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    axis.line = element_line(size = 1),  # 添加坐标轴线
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "",
       y = "NEP100") +
  scale_color_manual(values =c('#ce5c69','#5770a6'), name = NULL)  + # 设置箱式图的边框颜色
  scale_fill_manual(values = c('#ce5c69','#5770a6'),name = NULL)  + # 设置散点的填充色
  #scale_y_continuous(limits = c(-0.2, 0.3)) +
  geom_signif(comparisons = list(c('NEPC', 'ARPC')), test = "wilcox.test", map_signif_level = F,
              y_position = max(new.PN$NEP100))

#######spaital
expr.data <- Seurat::Read10X("G:\\importance\\undergraduated\\数据集整理\\single cell\\GSE230282\\GSE230282" )
seurat <- Seurat::CreateSeuratObject(counts = expr.data, project = 'GSE230282', assay = 'RNA')
seurat$slice <- 1
seurat$region <- 'anterior'

##标准化
seurat <- SCTransform(seurat)
img <- Seurat::Read10X_Image(image.dir = "G:/importance/undergraduated/数据集整理/single cell/GSE230282/img/")
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = seurat)]
seurat[['image']] <- img

seurat <- AddModuleScore(seurat,
                         features = list(NEP100.UP),
                         ctrl = 100,search = T,seed = 20000709,
                         name = 'NEP100.UP')
seurat <- AddModuleScore(seurat,
                         features = list(NEP100.DN),
                         ctrl = 100,search = T,seed = 20000709,
                         name = 'NEP100.DN')
new <- NEP100(seurat,type='single',layer='data',species = 'homo')

SpatialFeaturePlot(seurat,features = 'NEP100.UP1',pt.size=2,interactive = F)
SpatialFeaturePlot(seurat,features = 'NEP100.DN1',pt.size=2,interactive = F)
SpatialFeaturePlot(new,features = 'NEP100',pt.size=2,interactive = F)