#!/usr/bin/Rscript
############## 
# 1.a quick turtorial
if(F)
{
  install.packages("glmnet", repos = "https://cran.us.r-project.org")
  library(glmnet)
  
  data(QuickStartExample)
  x <- QuickStartExample$x
  y <- QuickStartExample$y
  
  fit <- glmnet(x, y)
  
  plot(fit)
  plot(fit,label = TRUE)
  print(fit)
  coef(fit, s = 0.1)
  
  set.seed(29)
  nx <- matrix(rnorm(5 * 20), 5, 20)
  predict(fit, newx = nx, s = c(0.1, 0.05))
  
  cvfit <- cv.glmnet(x, y)
  plot(cvfit)
  
  cvfit$lambda.min
  
  coef(cvfit, s = "lambda.min")
  
  predict(cvfit, newx = x[1:5,], s = "lambda.min")
  
  ########################
  
  wts <-  c(rep(1,50), rep(2,50))
  fit <- glmnet(x, y, alpha = 0.2, weights = wts, nlambda = 20)
  print(fit)
  
  fit$lambda
  
  fit <- glmnet(x, y)
  any(fit$lambda == 0.5)  # 0.5 not in original lambda sequence
  
  coef.apprx <- coef(fit, s = 0.5, exact = FALSE)
  coef.exact <- coef(fit, s = 0.5, exact = TRUE, x=x, y=y)
  cbind2(coef.exact[which(coef.exact != 0)],
         coef.apprx[which(coef.apprx != 0)])
  
  predict(fit, newx = x[1:5,], type = "response", s = 0.05)
  
  plot(fit, xvar = "lambda", label = TRUE)
  
  plot(fit, xvar = "dev", label = TRUE)
  
  cvfit <- cv.glmnet(x, y, type.measure = "mse", nfolds = 20)
  
  print(cvfit)
  
  cvfit$lambda.min
  
  predict(cvfit, newx = x[1:5,], s = "lambda.min")
  
  coef(cvfit, s = "lambda.min")
  
  foldid <- sample(1:10, size = length(y), replace = TRUE)
  cv1  <- cv.glmnet(x, y, foldid = foldid, alpha = 1)
  cv.5 <- cv.glmnet(x, y, foldid = foldid, alpha = 0.5)
  cv0  <- cv.glmnet(x, y, foldid = foldid, alpha = 0)
  
  par(mfrow = c(1,1))
  plot(cv1); plot(cv.5); plot(cv0)
  plot(log(cv1$lambda)   , cv1$cvm , pch = 19, col = "red",
       xlab = "log(Lambda)", ylab = cv1$name)
  points(log(cv.5$lambda), cv.5$cvm, pch = 19, col = "grey")
  points(log(cv0$lambda) , cv0$cvm , pch = 19, col = "blue")
  legend("topleft", legend = c("alpha= 1", "alpha= .5", "alpha 0"),
         pch = 19, col = c("red","grey","blue"))
  
}

############
# 2.data preparation
# options(timeout=120)
# options(repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# install.packages("reshape2")
# install.packages("dplyr")
# install.packages("tidyfst")
# install.packages("tidyr")
# install.packages("readxl")
# install.packages('caret')
# install.packages("scales")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")
# BiocManager::install("rtracklayer")
library(reshape2) 
library(dplyr)
library(tidyfst)
library(readxl)
library(tidyr)
library(biomaRt)
library(rtracklayer)
library(caret)
library(scales)

setwd('/home/user/data2/lit/project/spatialTrans/2018-science-zhuying')

rm(list=ls())

input <- read.table('macaque_mRNASeq/processed/nhp_development_RPKM_rmTechRep.txt')
# listEnsemblArchives()
my_mart <-useMart("ensembl",host = 'https://jul2016.archive.ensembl.org')
listDatasets(my_mart)
my_dataset <-  useDataset("hsapiens_gene_ensembl",
                          mart = my_mart)
# attributes <- listAttributes(my_dataset)
# dim(attributes)
# head(attributes)
gene_list <- getBM(attributes = c("ensembl_gene_id","gene_biotype","hgnc_symbol"),
                   mart = my_dataset)
unique(gene_list$gene_biotype)
pc_gene <- filter(gene_list,gene_biotype=='protein_coding')

for (i in c('HSB','RMB')) {
  RPKM <- dplyr::select(input, starts_with(i))
  # subset for test
  # RPKM <- RPKM[1:5,1:5]
  
  # subset for protein-coding genes
  ensembl_id <- gsub("\\..*$", "", rownames(RPKM))
  RPKM <- RPKM[c(ensembl_id %in% pc_gene$ensembl_gene_id),]
  
  # subset for brain_regions(NCX)
  RPKM <- data.frame(t(RPKM))
  sample_br <- data.frame(sample_br=rownames(RPKM))
  sample_br <- separate(sample_br,sample_br,c('sample',"brain_region"), sep = "\\.")
  NCX=c("MFC",
        "OFC",
        "DFC",
        "VFC",
        "M1C",
        "S1C",
        "IPC",
        "A1C",
        "STC", 
        "ITC",
        "V1C")
  RPKM <- RPKM[sample_br$brain_region %in% NCX,]
  sample_br <- sample_br[sample_br$brain_region %in% NCX,]
  
  # create y
  if (i=='HSB'){
    metadata <- read_excel('human_mRNASeq/metadata/mRNA-seq_Sample metadata.xlsx',sheet=1,na='NA')
    metadata <- as.data.frame(metadata)
    sample_days <- data.frame(sample=metadata[5:nrow(metadata),2],age=metadata[5:nrow(metadata),5])
  } else
  {
    metadata <- read_excel('macaque_mRNASeq/metadata/mRNA-seq_sample.metadata.xlsx',sheet=1,na='NA')
    metadata <- as.data.frame(metadata)
    sample_days <- data.frame(sample=metadata[3:(nrow(metadata)-2),2],
                                age=metadata[3:(nrow(metadata)-2),5])
  }
  
  sample_days <- sample_days[match(sample_br$sample,sample_days$sample),]
  sample_days <- mutate(sample_days,age=as.numeric(sample_days$age))

  assign(paste0("RPKM_",i),RPKM)
  assign(paste0("ensembl_id_",i),ensembl_id)
  assign(paste0("sample_br_",i),sample_br)
  assign(paste0("metadata_",i),metadata)
  assign(paste0("sample_days_",i),sample_days)
}

rm(RPKM,ensembl_id,sample_br,metadata,sample_days)

############
# 3.fit models

library(glmnet)

early_sample <- sample_days_HSB[sample_days_HSB$age>8*7 & sample_days_HSB$age <(20*365+38*7),][,'sample']
old_sample <- sample_days_HSB[sample_days_HSB$age>38*7 & sample_days_HSB$age <(60*365+38*7),][,'sample']
fit_sample_br <- data.frame(fit_sample_br=rownames(RPKM_HSB))
fit_sample_br <- separate(fit_sample_br,fit_sample_br,c('sample',"brain_region"), sep = "\\.")
early_RPKM_HSB <- as.matrix(RPKM_HSB[fit_sample_br$sample %in% early_sample,])
old_RPKM_HSB <- as.matrix(RPKM_HSB[fit_sample_br$sample %in% old_sample,])

fit_sample_br_early <- fit_sample_br[fit_sample_br$sample %in% early_sample,]
fit_sample_br_old <- fit_sample_br[fit_sample_br$sample %in% old_sample,]
days_early <- as.matrix(log(sample_days_HSB[match(fit_sample_br_early$sample,sample_days_HSB$sample),'age'],2))
days_old <- as.matrix(sqrt(sample_days_HSB[match(fit_sample_br_old$sample,sample_days_HSB$sample),'age']))

early_model <- cv.glmnet(early_RPKM_HSB, days_early, type.measure = "mse", nfolds = 10,alpha = 0.5)
late_model <- cv.glmnet(old_RPKM_HSB, days_old, type.measure = "mse", nfolds = 10,alpha = 0.5)
png(file="early_model.png")
plot(early_model)
dev.off()

png(file="late_model.png")
plot(late_model)
dev.off()
early_model$lambda.min
early_model$lambda.1se
late_model$lambda.min
late_model$lambda.1se

############
# 4.predict

RPKM_RMB <- as.matrix(RPKM_RMB)

days_predict <- data.frame(predict(early_model, newx = RPKM_RMB, s = "lambda.min"))

days_predict$sample_br <- row.names(days_predict)

days_predict <- separate(days_predict,sample_br,c('sample',"brain_region"), sep = "\\.")

days_predict %>% group_by(sample) %>% summarise(age.predicted=median(lambda.min)) -> days_predict

uniq <- unique(sample_days_RMB)

mutate(merge(days_predict,uniq,by="sample"),age.real=log(as.numeric(age),2)) %>% filter(age<2*365) -> res

mutate(res,age.predicted.2=2**age.predicted) ->res

MAE(res$age.predicted,res$age.real)

plot(res$age.real,res$age.predicted)

p <- ggplot(data = res, mapping = aes(x = age, y = age.predicted.2))

p + geom_point(color="chartreuse4") +
  scale_y_continuous(
    trans = log2_trans(),
    # breaks = trans_breaks("log2", function(x)2**x),
    breaks = c(64,113,128,256,512,1024),
    # labels = trans_format("log2", math_format(2^.x))
    labels = scales::number
  ) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = c(64,80,128,256,512),
    # labels = trans_format("log2", math_format(2^.x))
    labels = scales::number
  ) +
  geom_abline(slope = 1,intercept = 0,linetype = 2,size=0.5)+
  geom_vline(xintercept = 80,color="red",size=0.5,linetype=2)+
  geom_hline(yintercept = 113,color="red",size=0.5,linetype=2,show.legend = TRUE)








