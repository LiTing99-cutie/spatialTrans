##################
# handle with human

library(reshape2) 
library(dplyr)
library(tidyfst)
library(readxl)
library(tidyr)
library(biomaRt)
library(rtracklayer)
library(caret)

setwd('/home/user/data2/lit/project/spatialTrans/2018-science-zhuying')

rm(list=ls())

input <- read.table('macaque_mRNASeq/processed/nhp_development_RPKM_rmTechRep.txt')
# listEnsemblArchives()
my_mart <-useMart("ensembl",host = 'https://jul2016.archive.ensembl.org')
# listDatasets(my_mart)
my_dataset <-  useDataset("hsapiens_gene_ensembl",
                          mart = my_mart)
# attributes <- listAttributes(my_dataset)
# dim(attributes)
# head(attributes)
gene_list <- getBM(attributes = c("ensembl_gene_id","gene_biotype","hgnc_symbol"),
                   mart = my_dataset)
unique(gene_list$gene_biotype)
pc_gene <- filter(gene_list,gene_biotype=='protein_coding')

RPKM <- dplyr::select(input, starts_with('HSB'))

# subset for test
# RPKM <- RPKM[1:5,1:5]

# subset for protein-coding genes
RPKM <- mutate(RPKM,ensembl_gene_id=gsub("\\..*$", "", rownames(RPKM)))
RPKM <- merge(RPKM,pc_gene,by="ensembl_gene_id")
# remove all duplicates
RPKM %>% group_by(hgnc_symbol) %>% filter(n()==1) -> RPKM

# validation
rpkm <- read.table("/home/user/data/uplee/data/dataset/paper/2014-RNA-philipp_humanMacaque_PFC_RNASeq/GSE51264_rpkm_human.tsv")

# annotation file
gff <- readGFF("gencode.v16.annotation.gtf")
head(gff)

mapid <- gff[gff$type == "gene", c("gene_id", "gene_name","gene_type")]
head(mapid)
mapid_pc <- filter(mapid,gene_type=="protein_coding")

# subset for protein-coding genes
rpkm <- mutate(rpkm,ensembl_gene_id=rownames(rpkm))
rpkm <- merge(rpkm,mapid_pc,by.x="ensembl_gene_id",by.y="gene_id")
# remove all duplicates
rpkm %>% group_by(gene_name) %>% filter(n()==1) -> rpkm

rpkm <- merge(rpkm,data.frame(RPKM$hgnc_symbol),by.x="gene_name",by.y="RPKM.hgnc_symbol")
RPKM <- merge(RPKM,data.frame(rpkm$gene_name),by.x="hgnc_symbol",by.y="rpkm.gene_name")

rownames(rpkm) <- rpkm$gene_name
rpkm <- mutate(rpkm,gene_name=NULL,ensembl_gene_id=NULL,gene_type=NULL)

rownames(RPKM) <- RPKM$hgnc_symbol
RPKM <- mutate(RPKM,hgnc_symbol=NULL,ensembl_gene_id=NULL,gene_biotype=NULL)

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
metadata <- read_excel('human_mRNASeq/metadata/mRNA-seq_Sample metadata.xlsx',sheet=1,na='NA')
metadata <- as.data.frame(metadata)
sample_days <- data.frame(sample=metadata[5:nrow(metadata),2],age=metadata[5:nrow(metadata),5])
sample_days <- sample_days[match(sample_br$sample,sample_days$sample),]
sample_days <- mutate(sample_days,age=as.numeric(sample_days$age))


# fit models

library(glmnet)

early_sample <- sample_days[sample_days$age>8*7 & sample_days$age <(20*365+38*7),][,'sample']
old_sample <- sample_days[sample_days$age>38*7 & sample_days$age <(60*365+38*7),][,'sample']

sample_br <- data.frame(sample_br=rownames(RPKM))
sample_br <- separate(sample_br,sample_br,c('sample',"brain_region"), sep = "\\.")

early_RPKM <- as.matrix(RPKM[sample_br$sample %in% early_sample,])
old_RPKM <- as.matrix(RPKM[sample_br$sample %in% old_sample,])

sample_br_early <- sample_br[sample_br$sample %in% early_sample,]
sample_br_old <- sample_br[sample_br$sample %in% old_sample,]

days_early <- as.matrix(log(sample_days[match(sample_br_early$sample,sample_days$sample),'age'],2))
days_old <- as.matrix(sqrt(sample_days[match(sample_br_old$sample,sample_days$sample),'age']))

#####
# does 10-fold cross-validation based on mean squared error criterion
early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mse", nfolds = 10,alpha = 0.5)

if(F)
{
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
}

# predict

# meta <- read.table("/home/user/data/uplee/data/dataset/paper/2014-RNA-philipp_humanMacaque_PFC_RNASeq/philipp_rnaseq_meta.txt",
#                    fill = TRUE,header = TRUE)
# meta %>% mutate(age=Year*365+Day+38*7) %>% filter(age<12*365) %>% filter(grepl('human',Work_ID)) %>% 
#   dplyr::select("Work_ID","age") -> days_real


# rpkm <- data.frame(t(rpkm))
# rpkm <- rpkm[rownames(rpkm) %in% days_real$Work_ID,]
mae_all <- NULL
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
tmp <- merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID")
mae_all <- rbind(mae_all,data.frame(alpha=0.5,lambda.type='lambda.min',measure.type='mse',
                            mae=sum(abs(tmp$lambda.min-tmp$age.log))/nrow(tmp)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
tmp <- merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") 
mae_all <- rbind(mae_all,data.frame(alpha=0.5,lambda.type='lambda.1se',measure.type='mse',
                            mae=sum(abs(tmp$lambda.1se-tmp$age.log))/nrow(tmp)))

# does 10-fold cross-validation based on mean absolute error criterion
early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mae", nfolds = 10,alpha = 0.5)

days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
tmp <- merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") 
mae_all <- rbind(mae_all,data.frame(alpha=0.5,lambda.type='lambda.min',measure.type='mae',
                            mae=sum(abs(tmp$lambda.min-tmp$age.log))/nrow(tmp)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
tmp <- merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") 
mae_all <- rbind(mae_all,data.frame(alpha=0.5,lambda.type='lambda.1se',measure.type='mae',
                            mae=sum(abs(tmp$lambda.1se-tmp$age.log))/nrow(tmp)))

# alpha=0############## 
early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mse", nfolds = 10,alpha = 0)

days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=0,lambda.type='lambda.min',measure.type='mse',
                            mae=MAE(tmp$lambda.min,tmp$age.log)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=0,lambda.type='lambda.1se',measure.type='mse',
                            mae=MAE(tmp$lambda.1se,tmp$age.log)))

early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mae", nfolds = 10,alpha = 0)

days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=0,lambda.type='lambda.min',measure.type='mae',
                            mae=MAE(tmp$lambda.min,tmp$age.log)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=0,lambda.type='lambda.1se',measure.type='mae',
                            mae=MAE(tmp$lambda.1se,tmp$age.log)))

# alpha=1############## 
early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mse", nfolds = 10,alpha = 1)

days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=1,lambda.type='lambda.min',measure.type='mse',
                            mae=MAE(tmp$lambda.min,tmp$age.log)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mse',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=1,lambda.type='lambda.1se',measure.type='mse',
                            mae=MAE(tmp$lambda.1se,tmp$age.log)))

early_model <- cv.glmnet(early_RPKM, days_early, type.measure = "mae", nfolds = 10,alpha = 1)

days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.min"))
merge(mutate(days_predict,lambda.type='lambda.min',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=1,lambda.type='lambda.min',measure.type='mae',
                            mae=MAE(tmp$lambda.min,tmp$age.log)))
days_predict <- data.frame(predict(early_model, newx = as.matrix(rpkm), s = "lambda.1se"))
merge(mutate(days_predict,lambda.type='lambda.1se',measure.type='mae',sample=rownames(days_predict)),
      mutate(days_real,age.log=log(age,2)),by.x="sample",by.y="Work_ID") -> tmp
mae_all <- rbind(mae_all,data.frame(alpha=1,lambda.type='lambda.1se',measure.type='mae',
                            mae=MAE(tmp$lambda.1se,tmp$age.log)))



##################
# handle with macaque


