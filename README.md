# Base_barrier_cells

## About

The choroid plexus (ChP) is a highly understudied structure of the central nervous system (CNS). The structure hangs in the brain ventricles, is composed of an epithelial cell layer, which produces the cerebrospinal fluid (CSF) and forms the blood-CSF barrier. It encapsulates a stromal mix of fenestrated capillaries, fibroblasts and a broad range of immune cells. Here, we report that the ChP base region harbors unique fibroblasts that cluster together, are connected by tight junctions and seal the ChP stroma from brain and CSF, thereby forming ChP base barrier cells (ChP BBCs). ChP BBCs are derived from meningeal mesenchymal precursors, arrive early during embryonic development, are maintained throughout life and are conserved across species. Moreover, we provide transcriptional profiles and key markers to label ChP BBCs and observe a striking transcriptional similarity with meningeal arachnoid barrier cells (ABCs). Finally, we provide evidence that this fibroblast cluster functions as a barrier to control communication between CSF and the ChP stroma and between the latter and the brain parenchyma. Moreover, loss of barrier function was observed during an inflammatory insult. Altogether, we have identified a novel barrier that provides functional compartmentalization of ChP, brain and CSF.


## Overview scripts

Here's an overview of the various R scripts used in processing the scRNA-Seq data in the manuscript Verhaege et al.:
- 1.script_scRNAseq_LpsNegAggr.R: Seurat workflow for processing our 7-week-old ChP LV & 4V aggregate object in Figure 1
- 2.script_scRNAseq_GSM_2677817.R: Basic Seurat workflow for processing the ependymal sample 1 from Shah et al.
- 3.script_scRNAseq_GSM_2677818.R: Basic Seurat workflow for processing the ependymal sample 2 from Shah et al.
- 4.script_scRNAseq_GSM_2677819.R: Basic Seurat workflow for processing the ependymal sample 3 from Shah et al.
- 5.script_scRNAseq_Vanlandewijck.R: Basic Seurat workflow for processing the brain vascular cell data (GSE98816) from Vanlandewijck et al.
- 6.script_scRNAseq_Daan_and_Nina.R: Basic Seurat workflow for processing our 7- and 22-week-old ChP LV & 4V data
- 7.script_FB_datasets_Daan_loom.R: Data extraction and basic Seurat workflow for processing the CNS vascular and ependymal cell data from Zeisel et al. (2 loom files, l6_r3_vascular_cells and l6_r4_ependymal_cells, obtained from http://mousebrain.org/)
- 8.script_FB_datasets_Daan.R: Seurat and harmony workflow for merging our 7- and 22-week-old ChP LV & 4V data with various public datasets and as such creating the Fibroblast origin complete object in Figure 1
- 9.script_FB_datasets_Daan_clint.R: Further processing and exploration of the Fibroblast origin complete object in Figure 1
- 10.script_FB_datasets_Daan_subset2.R: Seurat and harmony workflow for merging the Fibroblasts from our 7- and 22-week-old ChP LV & 4V data with various public dataset Fibroblasts and as such creating the Fibroblast origin subset object in Figure 1
- 11.script_FB_datasets_Daan_subset2_clint.R: Further processing and exploration of the Fibroblast origin subset object in Figure 1
- 12.script_scRNAseq_Daan_and_Nina_v2.R: Basic Seurat workflow for processing our 7-, 22- and 82-week-old scRNA-Seq ChP LV & 4V data
- 13.script_scRNAseq_Lehtinen.R: Basic Seurat workflow for processing the embryonal scRNA-Seq ChP data from the Lehtinen lab (Dani et al.)
- 14.script_FB_datasets_Daan_Fig2_no_3V_CCA.R: Seurat and CCA workflow for merging the Fibroblasts from our 7-, 22- and 82-week-old scRNA-Seq ChP LV & 4V data with LV & 4V embryonal scRNA-Seq ChP data from the Lehtinen lab (Dani et al.) and as such creating the ChP Fibroblast age object in Figure 3. Further processing of the data also included in this script.
- 15.script_FB_datasets_human_CP.R: Exploration and attempted alternative processing of human snRNA-Seq data (Yang et al.) used in ChP Fibroblast species object in Figure 7
- 16.script_Harmony_snRNAseq_FBs_human_CP.R: Attempted Harmony analysis workflow snRNA-Seq human (Yang et al.) and snRNA-Seq mouse data from the Lehtinen lab (Dani et al.) for ChP Fibroblast species object in Figure 7
- 17.script_CCA_Integration_Urvb_and_snRNAseq_FBs_human_CP_Workflow1.R: Attempted CCA analysis workflow 1 snRNA-Seq human (Yang et al.) and scRNA-Seq mouse data for ChP Fibroblast species object in Figure 7
- 18.script_CCA_Integration_Urvb_and_snRNAseq_FBs_human_CP_Workflow2.R: Attempted CCA analysis workflow 2 snRNA-Seq human (Yang et al.) and scRNA-Seq mouse data for ChP Fibroblast species object in Figure 7
- 19.script_BBKNN_Integration_Urvb_and_snRNAseq_FBs_human_CP.R: BBKNN analysis snRNA-Seq human (Yang et al.) and scRNA-Seqmouse data for ChP Fibroblast species object in Figure 7

## Public data

1.	Shah, P. T. et al. Single-Cell Transcriptomics and Fate Mapping of Ependymal Cells Reveals an Absence of Neural Stem Cell Function. Cell 173, 1045-1057.e9 (2018).
2.	Vanlandewijck, M. et al. A molecular atlas of cell types and zonation in the brain vasculature. Nature 554, 475–480 (2018).
3.	DeSisto, J. et al. Single-Cell Transcriptomic Analyses of the Developing Meninges Reveal Meningeal Fibroblast Diversity and Function. Dev. Cell 54, 43-59.e4 (2020).
4.	Zeisel, A. et al. Molecular Architecture of the Mouse Nervous System. Cell 174, 999-1014.e22 (2018).
5.	Dani, N. et al. A cellular and spatial map of the choroid plexus across brain ventricles and ages. Cell 184, 3056-3074.e21 (2021).
6.	Yang, A. C. et al. Dysregulation of brain and choroid plexus cell types in severe COVID-19. Nat. 2021 5957868 595, 565–571 (2021).

## Citation

...
