library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratDisk)
library(data.table)
library(MetBrewer)

setwd("~/Downloads/")
allcase<-readRDS("snRNA_hc_aki_0627.rds")
ident<-allcase$orig.ident
dim(allcase)
counts <- GetAssayData(allcase, assay = "RNA")
counts2 <- counts[-(which(rownames(counts) %like% "MT-")),]
allcase<-CreateSeuratObject(counts2,min.cells = 1)
allcase$orig.ident<-ident

ident.list<-unique(allcase$orig.ident)
allcase$case.ident<-plyr::mapvalues(allcase$orig.ident,
                                    from=ident.list,
                                    to=c("3535","3535","3535","3535","3535","3504","3490","446","446",
                                         "460","460","461","462","462",ident.list[15:31]))
allcase$case.type.l1<-plyr::mapvalues(allcase$orig.ident,
                                      from=ident.list,
                                      to=c(rep("Healthy.Control",14),
                                           rep("AKI",17)))
allcase$case.type.l2<-plyr::mapvalues(allcase$orig.ident,
                                      from=ident.list,
                                      to=c(rep("Healthy.Control",14),
                                           rep("Non.COVID.AKI",13),
                                           rep("COVID.AKI",4)))

saveRDS(allcase,"snRNA_hc_aki_mtremoved_0705.rds")

seurat.list <- SplitObject(allcase, split.by = "case.ident")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
kidney.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")
allcase.combined <- IntegrateData(anchorset = kidney.anchors)


#PCA on integrated dataset
DefaultAssay(allcase.combined) <- "integrated"
allcase.combined <- ScaleData(allcase.combined, verbose = FALSE)
allcase.combined <- RunPCA(allcase.combined, npcs = 30, verbose = FALSE)
ElbowPlot(allcase.combined,ndims=30)

allcase.combined <- allcase.combined %>%FindNeighbors(dims=1:30,verbose=F)%>%FindClusters(resolution=0.5,verbose=F)%>%identity()
allcase.combined <- RunUMAP(allcase.combined, dims = 1:30,verbose=F)

DimPlot(allcase.combined, label=TRUE, pt.size = .1,label.size = 3,repel=T,split.by="case.type.l2")
feature=c("SLC5A12","SLC22A6","HAVCR1","LCN2","VCAM1","TACSTD2",
          "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
          'ACTA2',"DCN","MYH11", "EMCN","FLT1","CD163","IL7R")
DefaultAssay(allcase.combined) <- "RNA"
DotPlot(allcase.combined, features = feature,scale=T) +theme_classic()+ RotatedAxis()

#subset PT cells to recluster, because we need to put PT annotation back to the combined datasets. 
pt<-subset(allcase.combined,subset=integrated_snn_res.0.5==0|integrated_snn_res.0.5==3|integrated_snn_res.0.5==8|integrated_snn_res.0.5==12|integrated_snn_res.0.5==18)

#integrate and recluster PT subsets
seurat.list <- SplitObject(pt, split.by = "case.ident")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
pt.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")
pt.combined <- IntegrateData(anchorset = pt.anchors)

#PCA on integrated dataset
DefaultAssay(pt.combined) <- "integrated"
pt.combined <- ScaleData(pt.combined)
pt.combined <- RunPCA(pt.combined, npcs = 30, verbose = FALSE)
ElbowPlot(pt.combined,ndims=30)

DefaultAssay(pt.combined) <- "integrated"
i=10
pt.combined <- pt.combined %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.3,verbose=F)%>%identity()
pt.combined <- RunUMAP(pt.combined, dims = 1:i,verbose=F)
DimPlot(pt.combined,label=T)
DefaultAssay(pt.combined) <- "RNA"
DotPlot(pt.combined, features = c("CUBN","SLC5A12","SLC22A6","PRODH2","SLC7A13","HAVCR1","FTH1","FTL","SLC8A1","UMOD","SLC12A3",
                                  "LCN2","VCAM1","SOX4","CD24"),
        scale=T) +theme_classic()+ RotatedAxis()

#put PT annotation back to orignal annotation. 
DefaultAssay(pt.combined) <- "RNA"
pt.combined$integrated_snn_res.0.3<-plyr::mapvalues(pt.combined$integrated_snn_res.0.3,
                                                    from=c(3,1,4,0,5,2,7,6),
                                                    to=c("PT.S1S2","PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2",
                                                         "PT.Injured.Severe","PT.Injured.Severe","PT.Mixed"))

allcase.combined$integrated_snn_res.0.5.2<-as.character(allcase.combined$integrated_snn_res.0.5)
allcase.combined$integrated_snn_res.0.5.2[Cells(pt.combined)]<-as.character(pt.combined$integrated_snn_res.0.3)
DimPlot(allcase.combined,group.by = "integrated_snn_res.0.5.2",label=T,repel=T)

allcase.combined$celltype.orig<-plyr::mapvalues(allcase.combined$integrated_snn_res.0.5.2,
                                                from=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe","PT.Mixed",
                                                       11,16,2,1,14,6,
                                                       9,13,15,5,10,22,
                                                       21,7,19,
                                                       4,23,25,
                                                       17,20,24),
                                                to=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe","PT.Mixed",
                                                     "DTL","Degenerated","TAL.Injured","TAL.1","TAL.2","TAL.3",
                                                     "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                     "POD","Fib.1","Fib.2",
                                                     "EC","EC.GC","EC.AVR",
                                                     "Myeloid","TCell","BCell"))
allcase.combined$celltype.orig<-as.factor(allcase.combined$celltype.orig)
allcase.combined$celltype.orig<-factor(allcase.combined$celltype.orig,
                                       levels=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe","PT.Mixed",
                                                "DTL","Degenerated","TAL.Injured","TAL.1","TAL.2","TAL.3",
                                                "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                "POD","Fib.1","Fib.2",
                                                "EC","EC.GC","EC.AVR",
                                                "Myeloid","TCell","BCell"))

#subset TAL cells to recluster, because we need to put TAL annotation back to the combined datasets. 
tal<-subset(allcase.combined,subset=celltype.orig=="TAL.Injured"|celltype.orig=="TAL.1"|celltype.orig=="TAL.2"|celltype.orig=="TAL.3")
tal$case.ident2<-ifelse(tal$case.ident=="34_10050","32_2",
                        ifelse(tal$case.ident=="461","462",tal$case.ident))

seurat.list <- SplitObject(tal, split.by = "case.ident2")
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seurat.list)
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
kidney.anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features, reduction = "rpca")
tal.combined <- IntegrateData(anchorset = kidney.anchors)

#PCA on integrated TAL dataset
DefaultAssay(tal.combined) <- "integrated"
tal.combined <- ScaleData(tal.combined)
tal.combined <- RunPCA(tal.combined, npcs = 30, verbose = FALSE)
ElbowPlot(tal.combined,ndims=30)

DefaultAssay(tal.combined) <- "integrated"
i=10
tal.combined <- tal.combined %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.3,verbose=F)%>%identity()
tal.combined <- RunUMAP(tal.combined, dims = 1:i,verbose=F)
DimPlot(tal.combined,label=T)
DefaultAssay(tal.combined) <- "RNA"
DotPlot(tal.combined, features = c("UMOD","SLC12A1","NOS1","CLDN6","CLDN14","EGF","SLC12A3","SLC8A1","AQP2","FTH1","FTL",
                                   "HLA-A","HLA-C","HLA-E","CD74","HLA-DRA",
                                   "HAVCR1","LCN2","VCAM1","ROBO2","PROM1","SOX4","CD24","HNF1B"),
        scale=T) +theme_classic()+ RotatedAxis()
tal.combined$tal.cluster<-plyr::mapvalues(tal.combined$integrated_snn_res.0.3,
                                          from=c(5,0,4,7,9,
                                                 2,1,3,8,6),
                                          to=c("TAL.Injured","TAL.Injured","TAL.Injured","TAL.Injured","TAL.Mixed",
                                               "MD","TAL","TAL","TAL","TAL"))
DimPlot(tal.combined,label=T,group.by="tal.cluster")

#merge TAL cluster back to allcase
allcase.combined$celltype.orig<-as.character(allcase.combined$celltype.orig)
allcase.combined$celltype.orig[Cells(tal.combined)]<-as.character(tal.combined$tal.cluster)
allcase.combined$celltype.orig<-plyr::mapvalues(allcase.combined$celltype.orig,
                                                from=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe","PT.Mixed",
                                                       "DTL","Degenerated","TAL","TAL.Injured","TAL.Mixed","MD",
                                                       "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                       "POD","Fib.1","Fib.2",
                                                       "EC","EC.GC","EC.AVR",
                                                       "Myeloid","TCell","BCell"),
                                                to=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe","PT.Mixed",
                                                     "DTL","Degenerated","TAL","TAL.Injured","TAL.Mixed","MD",
                                                     "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                     "POD","Fib","vSMC.Per",
                                                     "EC","EC.Glom","EC.Fib",
                                                     "Myeloid","TCell","BCell"))

allcase.combined$celltype.orig<-factor(allcase.combined$celltype.orig,
                                       levels=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe","PT.Mixed",
                                                "DTL","Degenerated","TAL","TAL.Injured","TAL.Mixed","MD",
                                                "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                "POD","Fib","vSMC.Per",
                                                "EC","EC.Glom","EC.Fib",
                                                "Myeloid","TCell","BCell"))
allcase.combined<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed"&celltype.orig!="TAL.Mixed"&celltype.orig!="Degenerated")
