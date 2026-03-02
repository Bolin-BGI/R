# Utilize the function in seurat for data analysis
library(Seurat)
library(SoupX)
#library(scDblFinder)
#suppressPackageStartupMessages(library(GenomicRanges))
#suppressPackageStartupMessages(library(rtracklayer))
#函数list：
#(1):runSoupX_Droplets
#(2):runSoupX_NoDroplets
#(3):ExportData_10X_format
#(4):runDoubletFinder
#(5):runscDblFinder
#(6):runscDblFinder_MultiSamples
#(7):runAmulet

####################################### pre-processing #######################################
####################################### pre-processing #######################################
####################################### pre-processing #######################################

###1. remove Ambient RNAs by contamination
runSoupX_Droplets <- function(matrix_path=NULL, expression_matrix=NULL, droplet_matrix=NULL){
    # loading data
    if (!is.null(matrix_path)){
        toc = Seurat::Read10X(file.path(matrix_path, "04.Matrix", "FilterMatrix"), gene.column = 1)
        tod = Seurat::Read10X(file.path(matrix_path, "02.cDNAAnno", "RawMatrix"), gene.column = 1)
    }
    else{
        toc=expression_matrix
        tod=droplet_matrix
    }
    #
    genes=intersect(rownames(toc),rownames(tod))
    toc=toc[genes,]
    tod=tod[genes,]
    sc = SoupChannel(tod, toc)
    # get clusters info
    obj <- CreateSeuratObject(counts = toc)
    obj=Process_RNA(obj)
    soupx_groups = Idents(obj)
    #
    sc = setClusters(sc, soupx_groups)
    sc = autoEstCont(sc, doPlot=FALSE)
    out = adjustCounts(sc,roundToInt = TRUE)
    return (out)
}

runSoupX_NoDroplets <- function(matrix_path=NULL, expression_matrix=NULL){
    # loading data
    if (!is.null(matrix_path)){
        toc = Seurat::Read10X(matrix_path, gene.column = 1)
    }
    else{toc=expression_matrix}
    #
    scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
    # Calculate soup profile
    soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))
    scNoDrops = setSoupProfile(scNoDrops, soupProf)
    # get clusters info
    obj <- CreateSeuratObject(counts = toc)
    obj <- Process_RNA(obj)
    soupx_groups = Idents(obj)
    # Set cluster information in SoupChannel
    scNoDrops = setClusters(scNoDrops, soupx_groups)
    # Estimate contamination fraction
    scNoDrops  = autoEstCont(scNoDrops, doPlot=FALSE)
    # Infer corrected table of counts and rount to integer
    out = adjustCounts(scNoDrops, roundToInt = TRUE)
    return (out)
}

###2. convert seurat obj to 10X foramt
ExportData_10X_format <- function(data=NULL, out_dir=NULL, assay='RNA', if_compress=TRUE){
  #
  if (is.null(data)){stop('please supply the data')}
  #
  if(assay=='ATAC'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/peaks.bed'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
        R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
        R.utils::gzip(paste0(out_dir,'/peaks.bed'))
        R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  if(assay=='RNA'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/features.tsv'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
        R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
        R.utils::gzip(paste0(out_dir,'/features.tsv'))
        R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  print('Success to output data')
}


###3. identify doublets by DoubletFinder、scDblFinder、Amulet etc.

runDoubletFinder <- function(obj, dims, estDubRate=0.075, ncores=1){
  # Run DoubletFinder on a provided (preprocessed) Seurat object
  # Return the seurat object with the selected pANN parameter and the 
  # DoubletFinder doublet classifications

  ### pK Identification (parameter-sweep) ###
  # "pK ~ This defines the PC neighborhood size used to compute pANN (proportion of artificial nearest neighbors), 
  # expressed as a proportion of the merged real-artificial data. 
  # No default is set, as pK should be adjusted for each scRNA-seq dataset"

  sweep.res.list <- paramSweep_v3(obj, PCs=dims, sct=FALSE, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))

  # Get expected doublets (DF.classify will identify exactly this percent as doublets!)
  nExp_poi <- round(estDubRate * length(Cells(obj)))

  # DoubletFinder:
  obj <- doubletFinder_v3(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL

  return(obj)
}

# one sample/library
runscDblFinder <- function(obj, assay='ATAC'){
    #
    counts <- GetAssayData(obj[[assay]],slot = "counts")
    sce <- SingleCellExperiment(list(counts=counts))
    #
    sce <- scDblFinder(sce, clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
    dbl_score <- sce@colData[,c('scDblFinder.class','scDblFinder.score')]
    obj=AddMetaData(obj, metadata = as.data.frame(dbl_score[colnames(obj),c('scDblFinder.class', 'scDblFinder.score')]))
    #
    return (obj)
}

# multiple samples/libraries
runscDblFinder_MultiSamples <- function(obj, assay='ATAC', group='batch'){
    #
    counts <- GetAssayData(obj[[assay]],slot = "counts")
    sce <- SingleCellExperiment(list(counts=counts), colData=DataFrame(obj@meta.data[group]))
    #
    sce <- scDblFinder(sce, samples=group, clusters=TRUE, aggregateFeatures=TRUE, nfeatures=25, processing="normFeatures")
    dbl_score <- sce@colData[,c('scDblFinder.class','scDblFinder.score')]
    obj=AddMetaData(obj, metadata = as.data.frame(dbl_score[colnames(obj),c('scDblFinder.class', 'scDblFinder.score')]))    
    return (obj)
}

# for scATAC-seq data
runAmulet <- function(obj, fragment_path, species='human', save_path=NULL){
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(rtracklayer))
    #
    if (species=='human'){repeats <- blacklist_hg38}
    else {repeats <- blacklist_mm10}
    otherChroms <- GRanges(c("M","chrM","MT","X","Y","chrX","chrY"),IRanges(1L,width=10^8))
    toExclude <- suppressWarnings(c(repeats, otherChroms))
    #
    fragfile=fragment_path
    res <- amulet(fragfile, regionsToExclude=toExclude)
    # save Amulet results
    if (!is.null(save_path)){
        write.csv(res, save_path, quote=FALSE)
    }
    # add info into seurat obj
    obj=AddMetaData(obj, metadata = res[colnames(obj), c('p.value', 'q.value')])
    return (obj)
}

combineP_method <- function(obj){
    obj$scDblFinder.p <- 1-obj@meta.data[, "scDblFinder.score"]
    obj$combined <- apply(obj@meta.data[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
        x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
        suppressWarnings(aggregation::fisher(x))
    })
    #
    return (obj)
}

####################################### pre-processing #######################################
####################################### pre-processing #######################################
####################################### pre-processing #######################################

###1. run clusters using counts table
Process_RNA <- function(proj=proj, nfeatures=2000, pc.nums=50, dims.use=1:50, pca.name = 'pca', umap.name= 'umap', res = 0.8){
  DefaultAssay(proj) <- "RNA"
  #
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- ScaleData(proj)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = pca.name)
  proj <- FindNeighbors(proj, reduction = pca.name, dims = dims.use, verbose = FALSE)
  proj <- FindClusters(proj, resolution = res, algorithm = 3, verbose = FALSE)
  proj <- RunUMAP(proj, reduction = pca.name, dims = dims.use, reduction.name = umap.name)
  #
  return(proj)
}


###2. run clusters using counts table
Recluster_RNA <- function(proj=proj, nfeatures=3000, pc.nums=30, dims.use=1:30, pca.name = 'pca', umap.name= 'umap', res = 0.8){
  DefaultAssay(proj) <- names( proj@assays)[1]
  #
  #proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- ScaleData(proj)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = pca.name)
  proj <- FindNeighbors(proj, reduction = pca.name, dims = dims.use, verbose = FALSE)
  proj <- FindClusters(proj, resolution = res, algorithm = 3, verbose = FALSE)
  proj <- RunUMAP(proj, reduction = pca.name, dims = dims.use, reduction.name = umap.name)
  #
  return(proj)
}