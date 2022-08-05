#DEG for marker genes and FGSEA.
ptmarkers<-FindAllMarkers(pt.combined,min.pct = 0,logfc.threshold = 0)
write.csv(ptmarkers,"ptmarkers.07.14.csv")

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

for (i in c("PT.S3","PT.Redifferentiating","PT.Maladaptive",
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

#get DEG of PT.Maladaptive from AKI for TRIBE analysis
allcase.combined$celltype.orig<-as.character(allcase.combined$celltype.orig)
allcase.combined$celltype.orig2<-ifelse(allcase.combined$celltype.orig=="PT.S1S2"|allcase.combined$celltype.orig=="PT.S3",
                                        "PT.Healthy",allcase.combined$celltype.orig)
allcase.combined$celltype.orig2<-as.factor(allcase.combined$celltype.orig2)
Idents(allcase.combined)<-allcase.combined$celltype.orig2

allcase.combined2<-subset(allcase.combined,subset=celltype.orig!="PT.Mixed" & case.type.l1=="AKI")
pt.maladaptation.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Injured.2",min.pct = 0)
write.csv(pt.maladaptation.markers,"pt.maladaptation.markers.AKI.csv")
pt.healthy.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Healthy",min.pct = 0)
write.csv(pt.healthy.markers,"pt.health.markers.AKI.csv")

Idents(allcase.combined2)<-allcase.combined2$celltype.orig
pt.maladaptive.markers<-FindMarkers(allcase.combined2,ident.1 = "PT.Maladaptive",min.pct = 0)
write.csv(pt.maladaptive.markers,"pt.maladaptive.markers.AKI.csv")
