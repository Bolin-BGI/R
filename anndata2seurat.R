#' Convert AnnData to Seurat object
#'
#' This function converts an AnnData object (from Python) to a Seurat object (in R).
#'
#' @param inFile Input AnnData file path (character)  # AnnData 文件的完整路径
#' @param outFile Output Seurat object file path (character, NULL for no save) # 输出的 Seurat 对象保存路径。如果为 NULL，则不保存文件，仅返回 Seurat 对象。
#' @param main_layer Main data layer to convert ("counts", "data", or "scale.data") # 指定 AnnData 对象中要转换为主要数据层的名称，可选值为 "counts"（原始计数数据）、"data"（处理后的数据）或 "scale.data"（标准化后的数据）。
#' @param assay Name of the assay in Seurat object (default: "RNA") # 指定 Seurat 对象中 Assay 的名称。Assay 是 Seurat 中用于存储不同类型数据的结构，通常用于存储 RNA 数据。
#' @param use_seurat Whether to use Seurat's built-in ReadH5AD (logical)，是否调用 Seurat 包的 ReadH5AD 函数直接读取 AnnData 文件并转换为 Seurat 对象
#' @param lzf Whether to use LZF compression when reading AnnData (logical) # 是否使用 LZF 压缩格式读取 AnnData 文件。LZF 是一种快速压缩算法，适用于读取大型文件。
#' @param target_uns_keys List of specific keys to extract from AnnData's uns field # 指定需要从 AnnData 的 uns 字段中提取的键值对。uns 字段通常用于存储额外的元数据，如降维结果或聚类信息。
#'
#' @return Seurat object


# Example usage:
# source("anndata2seurat.R")
# seurat_obj <- anndata2seurat("input.h5ad", outFile = "output.rds", main_layer = "counts")


anndata2seurat <- function(inFile, outFile = NULL, main_layer = "counts", 
                           assay = "RNA", use_seurat = FALSE, lzf = FALSE, 
                           target_uns_keys = list()) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  
  main_layer <- match.arg(main_layer, c("counts", "data", "scale.data"))
  inFile <- path.expand(inFile)
  
  anndata <- reticulate::import("anndata", convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)
  
  if (use_seurat) {
    if (lzf) {
      tmpFile <- paste0(tools::file_path_sans_ext(inFile), ".decompressed.h5ad")
      ad <- anndata$read_h5ad(inFile)
      ad$write(tmpFile)
      tryCatch(
        {
          srt <- Seurat::ReadH5AD(tmpFile)
        },
        finally = {
          file.remove(tmpFile)
        }
      )
    } else {
      srt <- Seurat::ReadH5AD(inFile)
    }
  } else {
    ad <- anndata$read_h5ad(inFile)
    
    obs_df <- .obs2metadata(ad$obs)
    var_df <- .var2feature_metadata(ad$var)
    
    if (reticulate::py_to_r(sp$issparse(ad$X))) {
      X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(ad$X)))
    } else {
      X <- t(reticulate::py_to_r(ad$X))
    }
    colnames(X) <- rownames(obs_df)
    rownames(X) <- rownames(var_df)
    
    if (!is.null(reticulate::py_to_r(ad$raw))) {
      raw_var_df <- .var2feature_metadata(ad$raw$var)
      raw_X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(ad$raw$X)))
      colnames(raw_X) <- rownames(obs_df)
      rownames(raw_X) <- rownames(raw_var_df)
    } else {
      raw_var_df <- NULL
      raw_X <- NULL
    }
    
    if (main_layer == "scale.data" && !is.null(raw_X)) {
      assays <- list(Seurat::CreateAssayObject(data = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = "scale.data", new.data = X)
      message("X -> scale.data; raw.X -> data")
    } else if (main_layer == "data" && !is.null(raw_X)) {
      if (nrow(X) != nrow(raw_X)) {
        message("Raw layer was found with different number of genes than main layer, resizing X and raw.X to match dimensions")
        raw_X <- raw_X[rownames(raw_X) %in% rownames(X), , drop = FALSE]
        X <- X[rownames(raw_X), , drop = FALSE]
      }
      assays <- list(Seurat::CreateAssayObject(counts = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = "data", new.data = X)
      message("X -> data; raw.X -> counts")
    } else if (main_layer == "counts") {
      assays <- list(Seurat::CreateAssayObject(counts = X))
      message("X -> counts")
    } else {
      assays <- list(Seurat::CreateAssayObject(data = X))
      message("X -> data")
    }
    names(assays) <- assay
    Seurat::Key(assays[[assay]]) <- paste0(tolower(assay), "_")
    
    if (main_layer == "scale.data" && !is.null(raw_X)) {
      assays[[assay]]@meta.features <- raw_var_df
    } else {
      assays[[assay]]@meta.features <- var_df
    }
    
    project_name <- sub("\\.h5ad$", "", basename(inFile))
    srt <- new("Seurat", assays = assays, project.name = project_name, version = packageVersion("Seurat"))
    Seurat::DefaultAssay(srt) <- assay
    Seurat::Idents(srt) <- project_name
    
    srt@meta.data <- obs_df
    embed_names <- unlist(reticulate::py_to_r(ad$obsm_keys()))
    if (length(embed_names) > 0) {
      embeddings <- sapply(embed_names, function(x) as.matrix(reticulate::py_to_r(ad$obsm[x])), simplify = FALSE, USE.NAMES = TRUE)
      names(embeddings) <- embed_names
      for (name in embed_names) {
        rownames(embeddings[[name]]) <- colnames(assays[[assay]])
      }
      
      dim.reducs <- vector(mode = "list", length = length(embeddings))
      for (i in seq(length(embeddings))) {
        name <- embed_names[i]
        embed <- embeddings[[name]]
        key <- switch(name,
                      sub("_(.*)", "\\L\\1", sub("^X_", "", toupper(name)), perl = TRUE),
                      "X_pca" = "PC",
                      "X_tsne" = "tSNE",
                      "X_umap" = "UMAP"
        )
        colnames(embed) <- paste0(key, "_", seq(ncol(embed)))
        dim.reducs[[i]] <- Seurat::CreateDimReducObject(
          embeddings = embed,
          loadings = new("matrix"),
          assay = assay,
          stdev = numeric(0L),
          key = paste0(key, "_")
        )
      }
      names(dim.reducs) <- sub("X_", "", embed_names)
      
      for (name in names(dim.reducs)) {
        srt[[name]] <- dim.reducs[[name]]
      }
    }
  }
  
  srt@misc <- .uns2misc(ad, target_uns_keys = target_uns_keys)
  
  if (!is.null(outFile)) saveRDS(object = srt, file = outFile)
  
  return(srt)
}

#' Prepare cell metadata
#'
#' This function prepare cell metadata from AnnData.obs
#'
#' @param obs_pd Input AnnData.obs dataframe
#' @param assay Assay name, default "RNA" (str)
#'
#' @return AnnData object
#'
#' @import reticulate
.obs2metadata <- function(obs_pd, assay = "RNA") {
  obs_df <- .regularise_df(reticulate::py_to_r(obs_pd), drop_single_values = FALSE)
  attr(obs_df, "pandas.index") <- NULL
  colnames(obs_df) <- sub("n_counts", paste0("nCounts_", assay), colnames(obs_df))
  colnames(obs_df) <- sub("n_genes", paste0("nFeaturess_", assay), colnames(obs_df))
  return(obs_df)
}

#' Prepare feature metadata
#'
#' This function prepare feature metadata from AnnData.var
#'
#' @param var_pd Input AnnData.var dataframe
#'
#' @return AnnData object
#'
#' @import reticulate
.var2feature_metadata <- function(var_pd) {
  var_df <- .regularise_df(reticulate::py_to_r(var_pd), drop_single_values = FALSE)
  attr(var_df, "pandas.index") <- NULL
  colnames(var_df) <- sub("dispersions_norm", "mvp.dispersion.scaled", colnames(var_df))
  colnames(var_df) <- sub("dispersions", "mvp.dispersion", colnames(var_df))
  colnames(var_df) <- sub("means", "mvp.mean", colnames(var_df))
  colnames(var_df) <- sub("highly_variable", "highly.variable", colnames(var_df))
  return(var_df)
}


.uns2misc <- function(ad, target_uns_keys = list()) {
  uns_keys <- intersect(target_uns_keys, reticulate::py_to_r(ad$uns_keys()))
  misc <- sapply(uns_keys, function(x) reticulate::py_to_r(ad$uns[x]), simplify = FALSE, USE.NAMES = TRUE)
  return(misc)
}


#' Regularise dataframe
#'
#' This function checks if certain columns of a dataframe is of a single value
#' and drop them if required
#'
#' @param df Input data frame, usually cell metadata table (data.frame-like
#'   object)
#' @param drop_single_values Drop columns with only a single value (logical)
#'
#' @return Dataframe
.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[["name"]] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        paste("Dropping single category variables:"),
        paste(colnames(df)[k_singular], collapse = ", ")
      )
    }
    df <- df[, !k_singular, drop = F]
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
  }
  return(df)
}
