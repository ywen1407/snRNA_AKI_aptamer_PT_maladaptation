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
dim(allcase)
unique(allcase$orig.ident)

ident.list<-unique(allcase$orig.ident)
ident.list
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

allcase<-readRDS("snRNA_hc_aki_mtremoved_spurious__removed_0705.rds")
varlist<-c()
ident.list<-unique(allcase$orig.ident)
ident.list
for (i in 1:length(ident.list)){
  temp<-subset(allcase,subset=orig.ident==ident.list[i])
  temp <- NormalizeData(temp, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
  temp <- FindVariableFeatures(temp, selection.method = "vst", nfeatures =200,verbose = F)
  varlist<-c(varlist,temp@assays$RNA@var.features)
}
varlist<-unique(varlist)

allcase <- NormalizeData(allcase, normalization.method = "LogNormalize", scale.factor = 10000,verbose = F)
allcase@assays$RNA@var.features<-varlist
allcase<- ScaleData(allcase,verbose=F)
allcase<- RunPCA(allcase,verbose = F)
ElbowPlot(allcase,ndims=30)
i=15
allcase <- allcase %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.6,verbose=F)%>%identity()
allcase <- RunUMAP(allcase, dims = 1:i,verbose=F)
DimPlot(allcase, label=TRUE, pt.size = .1,label.size = 3,repel=T) #,split.by="case.type.l1"
table(allcase$case.ident,allcase$seurat_clusters)
ggsave(paste0(i,"_nuclei_dim_.png"),plot = last_plot(),scale = 1, width = 15, height = 5, dpi = 600)

feature=c("SLC5A12","SLC22A6","HAVCR1","LCN2","VCAM1","TACSTD2","AQP1",
          "UMOD","SLC12A1","NOS1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
          'ACTA2',"DCN","MYH11","EMCN","FLT1","CD163","IL7R","MKI67","TOP2A")
DotPlot(allcase, features = feature,scale=T) + RotatedAxis()
markers<-FindMarkers(allcase,ident.1 = 13)
markers2<-FindMarkers(allcase,ident.1 = 17)
markers3<-FindMarkers(allcase,ident.1 = 24)

allcase$seurat_clusters<-as.numeric(allcase$seurat_clusters)
allcase2<-subset(allcase,subset=seurat_clusters<24)
dim(allcase2)
saveRDS(allcase2,"snRNA_hc_aki_mtremoved_spurious__removed_0705.rds")

#try RPCA for batch correction
allcase<-readRDS("snRNA_hc_aki_mtremoved_spurious__removed_0705.rds")
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
save(kidney.anchors, file='rpca_kpmp_0706.RData')

load('rpca_kpmp__pt_0706.RData')
allcase.combined <- IntegrateData(anchorset = kidney.anchors)

saveRDS(allcase.combined,"rpca_kpmp_0706.rds")

#PCA on integrated dataset
DefaultAssay(allcase.combined) <- "integrated"
allcase.combined <- ScaleData(allcase.combined, verbose = FALSE)
allcase.combined <- RunPCA(allcase.combined, npcs = 30, verbose = FALSE)
ElbowPlot(allcase.combined,ndims=30)

#checking which PC should work well
#PC 16- 30 all work generally well.
for (i in c(14:30)){
  DefaultAssay(allcase.combined) <- "integrated"
  allcase.combined <- allcase.combined %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.5,verbose=F)%>%identity()
  allcase.combined <- RunUMAP(allcase.combined, dims = 1:i,verbose=F)
  DimPlot(allcase.combined, label=TRUE, pt.size = .1,label.size = 3,repel=T)
  ggsave(paste0(i,"_nuclei_dim.png"),plot = last_plot(),scale = 1, width = 6, height = 5, dpi = 600)
  DimPlot(allcase.combined, label=TRUE, pt.size = .1,label.size = 3,repel=T,split.by="case.type.l2")
  ggsave(paste0(i,"_nuclei_dim_split.png"),plot = last_plot(),scale = 1, width = 15, height = 5, dpi = 600)
  feature=c("SLC5A12","SLC22A6","HAVCR1","LCN2","VCAM1","TACSTD2",
            "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
            "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
            'ACTA2',"DCN","MYH11", "EMCN","FLT1","CD163","IL7R")
  DefaultAssay(allcase.combined) <- "RNA"
  DotPlot(allcase.combined, features = feature,scale=T) +theme_classic()+ RotatedAxis()
  ggsave(paste0(i,"_nuclei_dot.png"),plot = last_plot(),scale = 1, width = 10, height = 5, dpi = 600)
}

feature=c("CUBN","SLC5A12","SLC22A6","PRODH2","SLC7A13","HAVCR1","LCN2","VCAM1","TACSTD2","SLC44A5",
          "UMOD","SLC12A1","SLC12A3","AQP2","SCNN1G",
          "SLC26A7","SLC4A9","CRB2","CLDN1","NPHS2",
          'ACTA2',"DCN","MYH11","NTRK3","EMCN","FLT1","IGTA8","SERPINE2","PLVAP","CD163","IL7R","MS4A1")
DefaultAssay(allcase.combined) <- "RNA"
DotPlot(allcase.combined, features = feature,scale=T,group.by="celltype.orig") +theme_classic()+ RotatedAxis()

allcase.combined<-readRDS("rpca_kpmp_0706.rds")
allcase.combined<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed")
table(allcase.combined$case.ident,allcase.combined$case.type.l1)
DotPlot(allcase.combined,features = c("EGF","UMOD"),group.by = "celltype.orig")

pt<-subset(allcase.combined,subset=integrated_snn_res.0.5==0|integrated_snn_res.0.5==3|integrated_snn_res.0.5==8|integrated_snn_res.0.5==12|integrated_snn_res.0.5==18)
saveRDS(pt,"rpca_kpmp_pt_0706.rds")

#integrate and recluster PT subsets
pt.anchors<-get(load('rpca_kpmp_pt_0706.rdata'))
pt.combined <- IntegrateData(anchorset = pt.anchors)
saveRDS(pt.combined,"rpca_kpmp_pt_rpca_0706.rds")

#PCA on integrated dataset
pt.combined<-readRDS("rpca_kpmp_pt_0706_rpca.RDS")
DefaultAssay(pt.combined) <- "integrated"
pt.combined <- ScaleData(pt.combined)
pt.combined <- RunPCA(pt.combined, npcs = 30, verbose = FALSE)
ElbowPlot(pt.combined,ndims=30)
table(pt.combined$orig.ident)
DefaultAssay(pt.combined) <- "integrated"
i=10
pt.combined <- pt.combined %>%FindNeighbors(dims=1:i,verbose=F)%>%FindClusters(resolution=0.3,verbose=F)%>%identity()
pt.combined <- RunUMAP(pt.combined, dims = 1:i,verbose=F)
DimPlot(pt.combined,label=T)
DefaultAssay(pt.combined) <- "RNA"
DotPlot(pt.combined, features = c("CUBN","SLC5A12","SLC22A6","PRODH2","SLC7A13","HAVCR1","FTH1","FTL","SLC8A1","UMOD","SLC12A3",
                         "LCN2","VCAM1","SOX4","CD24"),
        scale=T) +theme_classic()+ RotatedAxis()

#remove mixed identify from PT for differential gene expression
DefaultAssay(pt.combined) <- "RNA"
pt.combined$integrated_snn_res.0.3<-plyr::mapvalues(pt.combined$integrated_snn_res.0.3,
                                           from=c(3,1,4,0,5,2,7,6),
                                           to=c("PT.S1S2","PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2",
                                                "PT.Injured.Severe","PT.Injured.Severe","PT.Mixed"))
Idents(pt.combined)<-pt.combined$integrated_snn_res.0.3
DimPlot(pt.combined,group.by = "integrated_snn_res.0.3",label=T)

allcase.combined$integrated_snn_res.0.5.2<-as.character(allcase.combined$integrated_snn_res.0.5)
allcase.combined$integrated_snn_res.0.5.2[Cells(pt.combined)]<-as.character(pt.combined$integrated_snn_res.0.3)
DimPlot(allcase.combined,group.by = "integrated_snn_res.0.5.2",label=T,repel=T)
table(allcase.combined$integrated_snn_res.0.5.2)

#get differentially expressed genes in PT cells
#PT.S1S2, PT.S3 are compared against injury PT as a whole cluster
#PT.Injury.1, PT.Injury.2, PT.Severe Injury are compared against healthy PT as a whole cluster
#comparing injured cell types to healthy cell types
pt.combined<-readRDS("rpca_kpmp_pt_rpca_0706.rds")
pt.combined<-subset(pt.combined,subset=integrated_snn_res.0.3!="PT.Mixed")

#DEG for marker genes and FGSEA.
Idents(pt.combined)<-pt.combined$integrated_snn_res.0.3
ptmarkers2<-FindAllMarkers(pt.combined,min.pct = 0,logfc.threshold = 0)
write.csv(ptmarkers2,"ptmarkers.07.14.csv")

DotPlot(pt.combined,features=("SOX4"),group.by="integrated_snn_res.0.3")
pt.combined$integrated_snn_res.0.3
#fgsea of PT subclusters
library(devtools)
library(fgsea)
library(ggrepel)

pathways <- gmtPathways("c5.go.bp.v7.4.symbols.gmt")
cluster="PT.S1S2"
temp=ptmarkers[ptmarkers$cluster==cluster,]
ranks <- temp$avg_log2FC
names(ranks) <- temp$gene
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = ranks,
                  minSize  = 5,
                  maxSize  = 500)
fgseaRes$cluster=cluster
fgseaRes=fgseaRes[order(pval),]

for (i in c("PT.S3","PT.Injured.1","PT.Injured.2",
            "PT.Injured.Severe")){
  cluster=i
  temp=ptmarkers[ptmarkers$cluster==cluster,]
  ranks <- temp$avg_log2FC
  names(ranks) <- temp$gene
  temp_res <- fgsea(pathways = pathways, stats    = ranks, minSize  = 5, maxSize  = 500)
  temp_res$cluster=cluster
  temp_res=temp_res[order(pval),]
  fgseaRes<-rbind(fgseaRes,temp_res)
}
fgseaRes <- apply(fgseaRes,2,as.character)
write.csv(fgseaRes,"fgsea_pt_0714.csv",row.names = F)

fgseaRes<-as.data.frame(fgseaRes)
fgseaRes2<-fgseaRes[fgseaRes$padj<0.05 & fgseaRes$NES>=0,]
write.csv(fgseaRes2,"fgsea_pt_sig_0714.csv",row.names = F)

#visualizing GO terms
fgseaRes2<-read.csv("fgsea_pt_0714.csv")
go.s1<-fgseaRes2[fgseaRes2$padj<0.05 & fgseaRes2$NES>=0 & fgseaRes2$cluster=="PT.S1S2","pathway"]
go.s3<-fgseaRes2[fgseaRes2$padj<0.05 & fgseaRes2$NES>=0 & fgseaRes2$cluster=="PT.S3","pathway"]
go.injured.1<-fgseaRes2[fgseaRes2$padj<0.05 & fgseaRes2$NES>=0 & fgseaRes2$cluster=="PT.Injured.1","pathway"]
go.injured.2<-fgseaRes2[fgseaRes2$padj<0.05 & fgseaRes2$NES>=0 & fgseaRes2$cluster=="PT.Injured.2","pathway"]
go.injured.severe<-fgseaRes2[fgseaRes2$padj<0.05 & fgseaRes2$NES>=0 & fgseaRes2$cluster=="PT.Injured.Severe","pathway"]
setdiff(go.injured.severe,c(go.injured.1,go.injured.2,go.s1,go.s3))

go.healthy.selected<-c("GOBP_ORGANIC_ACID_CATABOLIC_PROCESS","GOBP_ORGANIC_ACID_METABOLIC_PROCESS",
                       "GOBP_SMALL_MOLECULE_CATABOLIC_PROCESS","GOBP_REGULATION_OF_FATTY_ACID_OXIDATION",
                       "GOBP_DE_NOVO_NAD_BIOSYNTHETIC_PROCESS",
                       "GOBP_ORGANIC_ANION_TRANSPORT") # not much different betwen S1S2 and S3. 
go.injured.1.selected<-c("GOBP_KIDNEY_EPITHELIUM_DEVELOPMENT","GOBP_TISSUE_MORPHOGENESIS",
                         "GOBP_CELL_CELL_SIGNALING_BY_WNT","GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY",
                         "GOBP_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                         "GOBP_MYOTUBE_DIFFERENTIATION","GOBP_SMOOTH_MUSCLE_CELL_MIGRATION",
                         "GOBP_ACTOMYOSIN_STRUCTURE_ORGANIZATION",
                         "GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
                         "GOBP_MYELOID_CELL_DIFFERENTIATION","GOBP_LEUKOCYTE_DIFFERENTIATION",
                         "GOBP_NOTCH_SIGNALING_PATHWAY")
go.injured.2.selected<-c("GOBP_NEPHRON_EPITHELIUM_DEVELOPMENT",
                         
                         "GOBP_IMMUNE_SYSTEM_DEVELOPMENT","GOBP_LEUKOCYTE_MIGRATION","GOBP_MYELOID_CELL_DEVELOPMENT","GOBP_LYMPHOCYTE_ACTIVATION",
                         "GOBP_LEUKOCYTE_ADHESION_TO_VASCULAR_ENDOTHELIAL_CELL","GOBP_LEUKOCYTE_DIFFERENTIATION","GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
                         "GOBP_MYOTUBE_DIFFERENTIATION","GOBP_MYOBLAST_FUSION")
go.injured.severe.selected<-c("GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN","GOBP_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_PROTON_TRANSPORT",
                              
                              "GOBP_RESPONSE_TO_METAL_ION","GOBP_IRON_ION_TRANSPORT","GOBP_RESPONSE_TO_TOXIC_SUBSTANCE","GOBP_DETOXIFICATION_OF_INORGANIC_COMPOUND",
                              "GOBP_DNA_DAMAGE_RESPONSE_DETECTION_OF_DNA_DAMAGE",
                              "GOBP_HYDROGEN_PEROXIDE_CATABOLIC_PROCESS","GOBP_RESPONSE_TO_HYDROGEN_PEROXIDE",
                              "GOBP_RIBOSOMAL_LARGE_SUBUNIT_BIOGENESIS",
                              
                              "GOBP_REGULATION_OF_CELL_DEATH","GOBP_APOPTOTIC_SIGNALING_PATHWAY","GOBP_POSITIVE_REGULATION_OF_APOPTOTIC_SIGNALING_PATHWAY","GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY","GOBP_APOPTOTIC_MITOCHONDRIAL_CHANGES",
                              
                              "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY","GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
                              "GOBP_MYD88_DEPENDENT_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY","GOBP_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                              "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
                              "GOBP_NIK_NF_KAPPAB_SIGNALING",
                              "GOBP_HUMORAL_IMMUNE_RESPONSE",
                              
                              "GOBP_T_CELL_MEDIATED_CYTOTOXICITY","GOBP_REGULATION_OF_T_CELL_ACTIVATION","GOBP_POSITIVE_REGULATION_OF_T_CELL_PROLIFERATION",
                              "GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
                              "GOBP_AMYLOID_FIBRIL_FORMATION","GOBP_FIBROBLAST_PROLIFERATION")

goterm<-unique(c(go.healthy.selected,go.injured.1.selected,go.injured.2.selected,go.injured.severe.selected))
goterm<-unique(c(go.healthy.selected,
          "GOBP_KIDNEY_EPITHELIUM_DEVELOPMENT","GOBP_CELL_CELL_SIGNALING_BY_WNT","GOBP_NOTCH_SIGNALING_PATHWAY","GOBP_MYELOID_CELL_DIFFERENTIATION",
          "GOBP_MYELOID_CELL_DEVELOPMENT","GOBP_LYMPHOCYTE_ACTIVATION","GOBP_LEUKOCYTE_MIGRATION","GOBP_MYOTUBE_DIFFERENTIATION","GOBP_MYOBLAST_FUSION","GOBP_POSITIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION",
          
          "GOBP_RESPONSE_TO_METAL_ION","GOBP_RESPONSE_TO_TOXIC_SUBSTANCE","GOBP_RESPONSE_TO_HYDROGEN_PEROXIDE","GOBP_DNA_DAMAGE_RESPONSE_DETECTION_OF_DNA_DAMAGE","GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
          "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY","GOBP_T_CELL_MEDIATED_CYTOTOXICITY","GOBP_HUMORAL_IMMUNE_RESPONSE","GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
          "GOBP_NIK_NF_KAPPAB_SIGNALING","GOBP_NIK_NF_KAPPAB_SIGNALING","GOBP_AMYLOID_FIBRIL_FORMATION"
          ))

fsea2<-fgseaRes2[fgseaRes2$pathway %in% goterm,c("pathway","pval","padj","NES","cluster")]
fsea2$pathway<-as.factor(fsea2$pathway)
fsea2$pathway<-factor(fsea2$pathway,levels = goterm)
fsea2$cluster<-plyr::mapvalues(fsea2$cluster,from = c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2",
                                                      "PT.Injured.Severe"),
                               to=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive",
                                    "PT.Injured.Severe"))
fsea2$cluster<-as.factor(fsea2$cluster)
fsea2$cluster<-factor(fsea2$cluster,levels = c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive",
                                               "PT.Injured.Severe"))
fsea2$padj<-as.numeric(fsea2$padj)
fsea2$NES<-as.numeric(fsea2$NES)

fsea2 %>%
  filter(padj<0.05)%>%
  ggplot(aes(x=cluster,y=pathway,color=NES,size=-log10(padj)))+geom_point()+
  scale_y_discrete(limits=rev)+
  scale_color_gradientn(colors=rev(met.brewer("Cassatt1",type = "continuous")))+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",angle = 90,vjust = 0.3))+
  theme(axis.text.y = element_text(face="bold"))+
  xlab("Proximal Tubule Subcluster")+
  ylab("GO Biological Process")+
  scale_size(range=c(1,10),breaks = c(1,2,5,10))
ggsave("go.pt.subcluster.png",plot=last_plot(),scale=1,width = 8,height=8,dpi=600)

#get DEG of PT.Injured.2 from HC and AKI for ARIC and TRIBE analysis
allcase.combined<-readRDS("rpca_kpmp_0706.rds")
allcase.combined$celltype.orig<-as.character(allcase.combined$celltype.orig)
allcase.combined$celltype.orig2<-ifelse(allcase.combined$celltype.orig=="PT.S1S2"|allcase.combined$celltype.orig=="PT.S3",
                                        "PT.Healthy",allcase.combined$celltype.orig)
allcase.combined$celltype.orig2<-as.factor(allcase.combined$celltype.orig2)
Idents(allcase.combined)<-allcase.combined$celltype.orig2

allcase.combined2<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed" & case.type.l1=="Healthy.Control")
Idents(allcase.combined2)<-as.factor(allcase.combined2$celltype.orig2)
pt.maladaptation.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Injured.2",min.pct = 0)
write.csv(pt.maladaptation.markers,"pt.maladaptation.markers.HC.csv")
pt.healthy.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Healthy",min.pct = 0)
write.csv(pt.healthy.markers,"pt.health.markers.HC.csv")

allcase.combined2<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed" & case.type.l1=="AKI")
pt.maladaptation.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Injured.2",min.pct = 0)
write.csv(pt.maladaptation.markers,"pt.maladaptation.markers.AKI.csv")
pt.healthy.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Healthy",min.pct = 0)
write.csv(pt.healthy.markers,"pt.health.markers.AKI.csv")

Idents(allcase.combined2)<-allcase.combined2$celltype.orig
pt.degeneration.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Injured.Severe",min.pct = 0)
write.csv(pt.degeneration.markers,"pt.degeneration.markers.AKI.csv")

#save meta data and dimension reduction coordinates for RNA velocity
umap<-as.data.frame(pt.combined[["umap"]]@cell.embeddings)
write.csv(umap,"pt.cell_embeddings.csv")
write.csv(rownames(umap),"pt.cellID.csv")
write.csv(pt.combined@meta.data,"pt.metadata.csv")
write.csv(pt.combined$integrated_snn_res.0.3,"pt.clubsers.csv")

#save PT dataset as h5ad file for scenic
pt.combined$cluster.l1<-pt.combined$integrated_snn_res.0.3
pt.combined$cluster.l1<-as.character(pt.combined$cluster.l1)
pt.combined$case.ident<-as.character(pt.combined$case.ident)
pt.combined$case.type.l1<-as.character(pt.combined$case.type.l1)
pt.combined$case.type.l2<-as.character(pt.combined$case.type.l2)
SaveH5Seurat(pt.combined, filename = "pt.aki.exclude.mixed.h5Seurat",overwrite = T)
Convert("pt.aki.exclude.mixed.h5Seurat", dest = "h5ad")

##aptamer gene expression
pt.combined$integrated_snn_res.0.3<-factor(pt.combined$integrated_snn_res.0.3,levels=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2",
                                     "PT.Injured.Severe"))
Idents(pt.combined)<-pt.combined$integrated_snn_res.0.3
DotPlot(pt,features = c("FSTL3","NLGN4X","COL23A1","TGFB2",
                        "AFM","PLG","ENPP6","P4HA2","PDK2"),scale = T
        )+RotatedAxis()

DimPlot(allcase.combined,group.by="integrated_snn_res.0.5.2",label=T)
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
DimPlot(allcase.combined,group.by="celltype.orig",label=T,repel = T)

DotPlot(allcase.combined2,
        features = c("BSG","GNS","MDH1","LDHB","EIF4EBP2","PEBP1","CTSZ","HAVCR2","CALR","FABP3","HINT1","HSPA8",
                     "CLSTN3","STX3"),scale = T,group.by = "celltype.orig")+
  theme_classic()+
  RotatedAxis()
ggsave("aptamer_nuclei_dot_AKI.png",plot=last_plot(),scale=1,width = 10,height=5,dpi=600)

DotPlot(allcase.combined2,
        features = c("PLXDC2","P4HA2","GPR135","ARFGAP1","ENPP6","USH1C","ROR1","ERBB3",
                     "PROC","PLG","AFM","PTPRD"),scale = T,group.by = "celltype.orig")+
  theme_classic()+
  RotatedAxis()
ggsave("aptamer_nuclei_dot_healthy_AKI.png",plot=last_plot(),scale=1,width = 10,height=5,dpi=600)

DotPlot(allcase.combined,
        features = c("TGFB2","COL23A1","CD200","NLGN4X","P4HA2","ENPP6","PROC","PLG","AFM"),scale = T,group.by = "celltype.orig")+
  theme_classic()+
  RotatedAxis()
ggsave("aptamer_nuclei_dot_AKI.png",plot=last_plot(),scale=1,width = 10,height=5,dpi=600)
######################random plotting

#pyroptosis process:
#activate NLRP3, bind ASC -> activate caspase-1 -> cleave IL1b, IL18 -> release via GSDMD
#activate caspase 4,5,11 -> GSDMD pore -> release DAMPs
#activate caspase 3-> GSDME pore -> release DAMPs
DotPlot(pt.combined,
        features = c("NLRP3","ASC","CASP1","IL1B","IL18","GSDMD",
                     "CASP4","CASP5","CASP11","CASP3","GSDME"),scale = T,group.by = "integrated_snn_res.0.3")+
  theme_classic()+
  RotatedAxis()+
  scale_size(range=c(1,7),breaks = c(0,1,2,5,10))

DotPlot(pt.combined,
        features = c("SLC7A11","SLC3A2","GCLC","GCLM","GSS","GPX4","GSR","ACSL4","SOX9"),scale = T,group.by = "integrated_snn_res.0.3")+
  theme_classic()+
  RotatedAxis() #ferroptosis genes from Ide

#necroptosis
#PMID 30341423
#TNFR1, CD95 (TNFRSF6),TNFRSF10A,TNFRSF10B (TRAIL-R1, TRAIL-R2): death receptors
#TLR3, TLR4, ZBP1, DAI, DLM-1: triger activation of RIPk1
DotPlot(pt.combined,
        features = c("TNFRSF1A","FAS","TNFRSF10A","TNFRSF10B",
                     "TLR3","TLR4","ZBP1","RIPK1","RIPK3","MLKL"),scale = T,group.by = "integrated_snn_res.0.3")+
  theme_classic()+
  RotatedAxis()+
  scale_size(range=c(1,7)) #necroptosis genes from 

#apoptosis
DotPlot(pt.combined,
        features = c("CASP3","CASP7","CASP8","BAD","BAK1","BCAP31","BNIP3","CTSC","CTSH",
                     "DAPK3","EEF1E1","MTCH2","NUPR1","PPP2R1A"),scale = T,group.by = "integrated_snn_res.0.3")+
  theme_classic()+
  RotatedAxis()

DotPlot(pt.combined,
        features = c("DIABLO"),scale = T,group.by = "integrated_snn_res.0.3")+
  theme_classic()+
  RotatedAxis() #apoptosis genes from Balzer

##recluster TAL
allcase.combined<-readRDS("rpca_kpmp_0706.rds")
tal<-subset(allcase.combined,subset=celltype.orig=="TAL.Injured"|celltype.orig=="TAL.1"|celltype.orig=="TAL.2"|celltype.orig=="TAL.3"|celltype.orig=="DCT.Injured")
#31549 TAL cells.
#merge ID 30_10050 and 32_2 (two cases with small number of cells); 
#merge 461 (small number) to 462 prior to integration
table(tal$case.ident)
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
saveRDS(tal.combined,"rpca_kpmp_tal_rpca_0718.rds")

#PCA on integrated TAL dataset
tal.combined<-readRDS("rpca_kpmp_tal_rpca_0718.rds")
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
allcase.combined<-readRDS("rpca_kpmp_0706.rds")
unique(allcase.combined$celltype.orig)
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
DimPlot(allcase.combined,group.by="celltype.orig",label=T,repel = T)
saveRDS(allcase.combined,"rpca_kpmp_0706.rds")

#Plotting section for manuscript
#need to remove mixed identify first
allcase.combined<-readRDS("rpca_kpmp_0706.rds")
allcase.combined<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed"&celltype.orig!="TAL.Mixed"&celltype.orig!="Degenerated")

#DimPlot of all cells
DimPlot(allcase.combined,group.by="celltype.orig",label=T,repel = T)+
  scale_color_manual(values=met.brewer("Lakota",n=27))+
  ggtitle("Kidney Cell Type")
ggsave("kidney.dimplot.png",plot = last_plot(),scale = 1, width = 8, height = 8, dpi = 600)

DimPlot(allcase.combined,group.by="celltype.orig",label=F,repel = T)+
  scale_color_manual(values=met.brewer("Lakota",n=27))+
  ggtitle("Kidney Cell Type")

#Dotplot of all cells
DefaultAssay(allcase.combined) <- "RNA"
feature=c("CUBN","LRP2","SLC5A12","SLC22A7","SLC7A13","HAVCR1","VCAM1","FTH1","FTL","TACSTD2","SLC44A5",
          "UMOD","SLC12A1","NOS1","LCN2","SLC12A3","CALB1","AQP2","SCNN1G",
          "SLC4A1","SLC26A4","NPHS2",
          'ACTA2',"DCN","MYH11","NTRK3","EMCN","FLT1","HECW2","CD163","IL7R","MS4A1")
allcase.combined$celltype.orig.rev<-allcase.combined$celltype.orig
allcase.combined$celltype.orig.rev<-factor(allcase.combined$celltype.orig.rev,
                                           levels=rev(c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe","PT.Mixed",
                                                        "DTL","TAL","TAL.Injured","TAL.Mixed","MD",
                                                        "DCT","CNT","CD.PC","CD.PC.Injured","CD.ICA","CD.ICB",
                                                        "POD","Fib","vSMC.Per",
                                                        "EC","EC.Glom","EC.Fib",
                                                        "Myeloid","TCell","BCell")))
DotPlot(allcase.combined, features = feature,scale=T,group.by="celltype.orig.rev") +
  theme_classic()+ RotatedAxis()+
  scale_color_distiller(direction = 1)+
  xlab("Gene")+
  ylab("Major Kidney Epithelial, Vascular, Stromal and Immune Cell")
ggsave("kidney.dotplot.png",plot = last_plot(),scale = 1, width = 8, height = 6, dpi = 600)

#PT dimplot
pt.combined<-readRDS("rpca_kpmp_pt_rpca_0706.rds")
pt.combined<-subset(pt.combined,subset=integrated_snn_res.0.3!="PT.Mixed")
pt.combined$celltype.orig<-plyr::mapvalues(pt.combined$integrated_snn_res.0.3,
                                                from=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe"),
                                                to=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe"))
pt.combined$celltype.orig<-factor(pt.combined$celltype.orig,
                                  levels=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe"))
DimPlot(pt.combined,group.by="celltype.orig",label=T,repel = T)+
  scale_color_manual(values=met.brewer("Cross",n=6))+
  ggtitle("Proximal Tubule Subcluster")
ggsave("pt.dimplot.png",plot = last_plot(),scale = 1, width = 6, height = 6, dpi = 600)

#PT dotplot
pt.combined$celltype.orig.rev<-pt.combined$celltype.orig
pt.combined$celltype.orig.rev<-factor(pt.combined$celltype.orig.rev,
                                  levels=rev(c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe")))

DotPlot(pt.combined, features = c("CUBN","LRP2","SLC5A12","SLC22A6","SLC22A7","SLC7A13","SPP1","FTH1","FTL","HLA-A","HLA-C","HLA-E","CD74","HLA-DRA",
                                  "HAVCR1","VCAM1","DCDC2","ROBO2","PROM1","SOX4","CD24"),group.by="celltype.orig.rev",
        scale=T) +theme_classic()+ RotatedAxis()+
  scale_color_distiller(direction = 1)+
  xlab("Gene")+
  ylab("Proximal Tubule Subcluster")
ggsave("pt.dotplot.png",plot = last_plot(),scale = 1, width = 8, height = 5, dpi = 600)

#ferroptosis
#PMID: 35803244
#protectors: 
#SLC7A11, SLC3A2 (importer of cysteine), GCLC (generate reduced glutathione)
#GPX4 (use reduced glutathione to reduce reactive PUFA-PL-OOH)
#iPLA2B (PLA2G6): similar function as GPX4
#prom2: pump reduced iron back to cytoplasm
#ACLS3, IL4i1, FSP1(gene name S100A4), DHODH, GCH1: suppress ALOX, PEBP1, POR

#damager:
#ALOXs, PEBP1, POR: convert PUFA-PL to toxic form under reduced iron.
#ALOX: ALOX12, ALOX12B, ALOXE3
#ACSL4, LPCAT3, ACSL1: convert Acetyl-CoA to PUFA-PL
#PKCbeta2 further activates ACSL4

#mixed results: some protectors of ferroptosis is activated in severe injury PT (GPX4) but lost in injured 1.2 
#certain damagers are activated even at healthy states, such as POR, ACSL1, LPCAT3. 
#but ACSL4 is activated in injured states 2. 
#different therapeutic targets in acute phase and aki-to-ckd transition (35210424)
DotPlot(pt.combined,
        features = c("SLC7A11","SLC3A2","GCLC","GCLM","GSS","GPX4","ACSL4"),scale = T,group.by = "celltype.orig.rev")+
  theme_classic()+
  RotatedAxis()+
  scale_size(range=c(1,7),breaks = c(0,5,10,25,50))+
  #scale_color_distiller(direction = 1)+
  scale_color_gradientn(colors=rev(met.brewer("Benedictus",type="continuous")))+
  xlab("Genes Involved in Ferroptosis")+
  ylab("Proximal Tubule Subcluster")
ggsave("pt.ferroptosis.png",plot = last_plot(),scale = 1, width = 8, height = 4, dpi = 600)


#necroptosis
DotPlot(pt.combined,
        features = c("TNFRSF1A","FAS","TNFRSF10A","TNFRSF10B",
                     "TLR3","TLR4","ZBP1","RIPK1","RIPK3","MLKL"),scale = T,group.by = "celltype.orig.rev")+
  theme_classic()+
  RotatedAxis()+
  #scale_color_distiller(direction = 1)+
  scale_color_gradientn(colors=rev(met.brewer("Benedictus",type="continuous")))+
  xlab("Genes Involved in Necroptosis")+
  ylab("Proximal Tubule Subcluster")+
  scale_size(range=c(1,7))
ggsave("pt.necroptosis.png",plot = last_plot(),scale = 1, width = 8, height = 4, dpi = 600)

#plotting PT composition
pt.combined$integrated_snn_res.0.3<-as.character(pt.combined$integrated_snn_res.0.3)
pt_table<-as.data.frame(cbind(pt.combined$integrated_snn_res.0.3,pt.combined$case.ident))
colnames(pt_table)<-c("cluster","case.id")
pt_table$count<-1
pt_table$cluster<-plyr::mapvalues(pt_table$cluster,
                                  from=c("PT.S1S2","PT.S3","PT.Injured.1","PT.Injured.2","PT.Injured.Severe"),
                                  to=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe"))
pt_table$cluster<-as.factor(pt_table$cluster)
pt_table$cluster<-factor(pt_table$cluster,levels=c("PT.S1S2","PT.S3","PT.Redifferentiating","PT.Maladaptive","PT.Injured.Severe"))
pt_table$case.id<-plyr::mapvalues(pt_table$case.id,
                                  from=c("30_10018","30_10034","30_10044","32_2",
                                         "32_10003","32_10034","32_10205","33_10005",
                                         "33_10006","34_10050","34_10184","34_10187",
                                         "34_10209",
                                         "cov_02_0002","cov_02_0010","cov_02_0094","cov_02_0096",
                                         "3490","3535","3504","446",
                                         "460","461",
                                         "462"),
                                  to=c("non.COVID.AKI.1","non.COVID.AKI.2","non.COVID.AKI.3","non.COVID.AKI.4",
                                       "non.COVID.AKI.5","non.COVID.AKI.6","non.COVID.AKI.7","non.COVID.AKI.8",
                                       "non.COVID.AKI.9","non.COVID.AKI.10","non.COVID.AKI.11","non.COVID.AKI.12",
                                       "non.COVID.AKI.13",
                                       "COVID.AKI.1","COVID.AKI.2","COVID.AKI.3","COVID.AKI.4",
                                       "Reference.1","Reference.2","Reference.3","Reference.4",
                                       "Reference.5","Reference.6",
                                       "Reference.7"))
pt_table$case.id<-as.factor(pt_table$case.id)
pt_table$case.id<-factor(pt_table$case.id,
                         levels=c("non.COVID.AKI.1","non.COVID.AKI.2","non.COVID.AKI.3","non.COVID.AKI.4",
                                  "non.COVID.AKI.5","non.COVID.AKI.6","non.COVID.AKI.7","non.COVID.AKI.8",
                                  "non.COVID.AKI.9","non.COVID.AKI.10","non.COVID.AKI.11","non.COVID.AKI.12",
                                  "non.COVID.AKI.13",
                                  "COVID.AKI.1","COVID.AKI.2","COVID.AKI.3","COVID.AKI.4",
                                  "Reference.1","Reference.2","Reference.3","Reference.4",
                                  "Reference.5","Reference.6",
                                  "Reference.7"))
ggplot(pt_table, aes(x = case.id, y=count,fill = cluster)) +
  geom_bar(position = "fill",stat="identity")+
  xlab("Participant ID")+
  ylab("Proportion")+
  theme_classic()+
  theme(axis.text.x = element_text(face="bold",angle = 90))+
  scale_fill_manual(values = met.brewer("Hiroshige",5))
ggsave("pt.composition.png",plot=last_plot(),scale=1,width = 7,height=5,dpi=600)

#aptamer dotplot
allcase.combined<-subset(allcase.combined,subset=case.type.l1=="AKI")
DotPlot(allcase.combined, features = c("NLGN4X","COL23A1","TGFB2","CD200","ENPP6","PLG","PROC",'P4HA2',"AFM"),
        scale=T,group.by="celltype.orig.rev") +
  theme_classic()+ RotatedAxis()+
  scale_color_distiller(direction = 1)+
  xlab("Gene for Biomarkers of PT Maladaptation and PT at Healthy States")+
  ylab("Major Kidney Epithelial, Vascular, Stromal and Immune Cell")
ggsave("pt.aptamer.dotplot.png",plot = last_plot(),scale = 1, width = 6, height = 6, dpi = 600)

#maladaptation candidate aptamer dotplot
DotPlot(allcase.combined, features = c("RRAS2","MPP6","CRIM1","EPHA7","BCL2","NLGN4X","TGFB2","COL23A1","CD200",'PDE7A',"PDE1A","ADAMTS3"),
        scale=T,group.by="celltype.orig.rev") +
  theme_classic()+ RotatedAxis()+
  scale_color_distiller(direction = 1)+
  xlab("Genes for Candidate Biomarkers of PT Maladaptation")+
  ylab("Major Kidney Epithelial, Vascular, Stromal and Immune Cell")
ggsave("pt.maladaptation.candidate.aptamer.dotplot.png",plot = last_plot(),scale = 1, width = 8, height = 6, dpi = 600)

#health candidate aptamer dotplot
DotPlot(allcase.combined, features = c("PLXDC2","GPR135","ARFGAP1","USH1C","ROR1","ERBB3","PTPRD","P4HA2",
                                       "ENPP6",'PLG',"PROC","AFM"),
        scale=T,group.by="celltype.orig.rev") +
  theme_classic()+ RotatedAxis()+
  scale_color_distiller(direction = 1)+
  xlab("Genes for Candidate Biomarkers of PT at Healthy States")+
  ylab("Major Kidney Epithelial, Vascular, Stromal and Immune Cell")
ggsave("pt.health.candidate.aptamer.dotplot.png",plot = last_plot(),scale = 1, width = 8, height = 6, dpi = 600)
