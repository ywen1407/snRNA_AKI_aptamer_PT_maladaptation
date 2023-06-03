# snRNA_AKI_aptamer_PT_maladaptation

This is the code for manuscript titled "Analysis of the Human Kidney Transcriptome and Plasma Proteome Identifies Novel Biomarkers of Proximal Tubule Maladaptation to Injury".

Author lists:

Authors: Yumeng Wen1, Emily Su2, Leyuan Xu3, Steven Menez1, Dennis Moledina3, Paul M. Palevsky4, Lloyd Cantley3, Patrick Cahan2, Chirag R. Parikh1*, for the Kidney Precision Medicine Project (KPMP) and Translational Investigation of Biomarker Endpoint of Acute Kidney Injury (TRIBE-AKI) Consortia

Affiliations:

Department of Medicine/Division of Nephrology, Johns Hopkins University School of Medicine, Baltimore, MD, USA.
Department of Biomedical Engineering, Johns Hopkins University School of Medicine, Baltimore, MD, USA.
Department of Medicine/Section of Nephrology, Yale School of Medicine, New Haven, CT, USA.
Renal-Electrolyte Division, University of Pittsburgh School of Medicine, Pittsburgh, PA, USA.
*Corresponding author. Email: chirag.parikh@jhmi.edu; Telephone: 410-955-5268

Datasets can be obtained by contacting the KPMP consortium; If you have any questions about the code, please contact first author at ywen14@jhmi.edu.

The general work flow is:
1. sequence alignment and library generation: cellranger v7
2. remove ambient RNA using cellbender
3. remove doublet using DoubletDetection
4. Preprocessing of combined datasets: including RPCA for data integration. remove clusters that express marker genes of more than 2 cell types
5. Select PT cluster: redo RPCA for data integration, remove one cluster that is potential doublet. 
6. Select TAL cluster: repeat step 5. 
7. DEG of PT subcluter in AKI samples only, then perform GSEA using FGSEA package.
8. Select PT subcluster using Scanpy, create h5ad dataset, use Epoch to reconstruct GRN topology, and use SCENIC to reconstruct cis-regulatory network and to calculate cellular enrichment. 
9. Integrate DEG of PT subcluster in AKI samples with proteomics findings in TRIBE-AKI. 
