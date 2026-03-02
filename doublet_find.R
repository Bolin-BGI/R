# /hwfssz5/ST_SUPERCELLS/P21Z10200N0171/USER/wangtao/script_hub/doublet_find.r
library(dplyr)
library(Seurat)
library(DoubletFinder)

Find_doublet <- function(data,ratio = 0.075){
        sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        nExp_poi <- round(as.numeric(ratio)*ncol(data))
        p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
        data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
        colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
        data
}