library(Seurat)
library(SeuratDisk)
library(dplyr)
library(epoch)

rm(list=ls())
gc()
memory.limit(size = 9999999999)

setwd("~/Downloads/")
pt.combined<-readRDS("rpca_kpmp_pt_rpca_0706.rds")
pt.combined<-subset(pt.combined,subset=case.type.l1=="AKI" & integrated_snn_res.0.3!="PT.Mixed")
hvg<-read.csv("hvg_4757.csv") # this is created from scanpy, using dispersion of 0.3. 
varlist<-hvg[hvg$highly_variable=="True","X_index"]
pt.combined@assays$RNA@var.features<-varlist

TF<-utils_loadObject("hsTFs.rda")
expDat<-as.matrix(pt.combined[["RNA"]]@data)
sampTab<-pt.combined@meta.data
expDat2<-expDat[rowSums(expDat)>0,]
dim(expDat2)
dgenes<-varlist
grnDF <- reconstructGRN(expDat2, TF, dgenes, method="pearson", zThresh=3)
grnDF<-grnDF[,c(2,1,3,4)]
grnDF<-grnDF[,c(1,2,3)]
colnames(grnDF)<-c("TF","target","importance")
write.csv(grnDF,"grnDF_epoch_pearson_pt_0715_4757hvg.csv",row.names = F)
