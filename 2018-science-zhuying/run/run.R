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
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library(reshape2) 
library(dplyr)
library(tidyfst)
library(readxl)
library(tidyr)
library(biomaRt)

setwd('/home/user/data2/lit/project/spatialTrans/2018-science-zhuying')

input <- read.table('macaque_mRNASeq/processed/nhp_development_RPKM_rmTechRep.txt')

RPKM <- input

RPKM <- dplyr::select(RPKM, starts_with("HSB")) 

# subset for test
# RPKM <- RPKM[1:5,1:5]

#subset for protein-coding genes
my_mart <-useMart("ensembl")
my_dataset <-  useDataset("hsapiens_gene_ensembl",
                          mart = my_mart)
gene_list <- getBM(attributes = c("ensembl_gene_id","gene_biotype"),
                   mart = my_dataset)
unique(gene_list$gene_biotype)
pc_gene <- filter(gene_list,gene_biotype=='protein_coding')

ensembl_id <- gsub("\\..*$", "", rownames(RPKM))
RPKM <- RPKM[c(ensembl_id %in% pc_gene$ensembl_gene_id),]

#subset for brain_regions
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

#select for NCX
RPKM <- RPKM[sample_br$brain_region %in% NCX,]

sample_br <- sample_br[sample_br$brain_region %in% NCX,]

metadata <- read_excel('human_mRNASeq/metadata/mRNA-seq_Sample metadata.xlsx',sheet=1,na='NA')
metadata <- as.data.frame(metadata)
sample_days <- data.frame(sample=metadata[5:nrow(metadata),2],age=metadata[5:nrow(metadata),5])

days <- matrix(log(as.numeric(sample_days[match(sample_br$sample,sample_days$sample),2]),2))

############
# 3.fit a model
RPKM <- as.matrix(RPKM)
library(glmnet)
cvfit <- cv.glmnet(RPKM, days, type.measure = "mse", nfolds = 10,alpha = 0.5)
print(cvfit)
png(file="cvfit.png")
plot(cvfit)
dev.off()

save(cvfit,file="cvfit.Rdata")

cvfit$lambda.min
cvfit$lambda.1se

############
# 4. predict
RPKM_r <- input

RPKM_r <- dplyr::select(RPKM_r, starts_with("RMB")) 

# subset for test
# RPKM_r <- RPKM_r[1:5,1:5]

#subset for protein-coding genes
ensembl_id_r <- gsub("\\..*$", "", rownames(RPKM_r))
RPKM_r <- RPKM_r[c(ensembl_id_r %in% pc_gene$ensembl_gene_id),]

#subset for brain_regions

RPKM_r <- data.frame(t(RPKM_r))

sample_br_r <- data.frame(sample_br_r=rownames(RPKM_r))

sample_br_r <- separate(sample_br_r,sample_br_r,c('sample',"brain_region"), sep = "\\.")

#select for NCX
RPKM_r <- RPKM_r[sample_br_r$brain_region %in% NCX,]

sample_br_r <- sample_br_r[sample_br_r$brain_region %in% NCX,]

metadata_r <- read_excel('macaque_mRNASeq/metadata/mRNA-seq_sample.metadata.xlsx',sheet=1,na='NA')
metadata_r <- as.data.frame(metadata_r)
sample_days_r <- data.frame(sample=metadata_r[3:(nrow(metadata_r)-2),2],
                            age=metadata_r[3:(nrow(metadata_r)-2),5])

days_r <- matrix(log(as.numeric(sample_days_r[match(sample_br_r$sample,sample_days_r$sample),2]),2))

RPKM_r <- as.matrix(RPKM_r)

days_r_predict <- data.frame(predict(cvfit, newx = RPKM_r, s = "lambda.min"))

days_r_predict$sample_br <- row.names(days_r_predict)

days_r_predict <- separate(days_r_predict,sample_br,c('sample',"brain_region"), sep = "\\.")

days_r_predict %>% group_by(sample) %>% summarise(age_predicted_log=median(lambda.min)) -> days_r_predict

days_r_predict %>%  mutate(age_predicted=2**age_predicted_log) -> days_r_predict

merge(days_r_predict,sample_days_r,by="sample") -> res_r

res_r %>% mutate(age_log=log(as.numeric(age),2)) -> res_r

plot(res_r$age_log,res_r$age_predicted_log)

# nohup Rscript run.R &>run.log &