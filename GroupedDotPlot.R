#' Seurat对象的按组分面气泡图 (Grouped DotPlot)
#'
#' 此函数利用ggplot2从Seurat对象生成分面气泡图（DotPlot）。
#' 允许用户输入一个包含基因集的命名列表，以自动对基因进行分组和分面，
#' 生成美观且达到出版级别的可视化图表。
#'
#' @param obj 一个Seurat对象。
#' @param gene_groups 一个命名列表（named list），其中名称为分组/分面标签，值为包含基因名称的向量。
#' @param group.by Seurat对象中用于映射到Y轴的元数据列名（例如：细胞类型 celltype）。
#' @param colors 一个包含3种颜色的向量，用于渐变色（低表达、中表达、高表达）。
#' @param size_range 长度为2的数值向量，表示气泡点大小的最小值和最大值。
#' @param base_text_size 图片主题的基础字体大小。
#' @param save_prefix 可选。如果提供该参数，则会自动将图片保存为PNG和PDF格式，例如："my_plot"（不需要手动加后缀）。
#' @param width 保存图表的宽度。
#' @param height 保存图表的高度。
#' 
#' @return 返回一个ggplot对象。
#' @export
#'
#' @import Seurat
#' @import ggplot2

GroupedDotPlot <- function(obj, 
                           gene_groups, 
                           group.by = "seurat_clusters",
                           colors = c("lightgrey", "#ffc7c7", "red"),
                           size_range = c(0.1, 5),
                           base_text_size = 12,
                           save_prefix = NULL,
                           width = 9, 
                           height = 8) {
  
  # 1. 将命名列表转换为映射数据框
  # 这一步是为了保留用户在列表中定义的分组和基因顺序
  group_names <- names(gene_groups)
  use_markers_df <- data.frame(
    features.plot = unlist(gene_groups, use.names = FALSE),
    type = rep(group_names, lengths(gene_groups))
  )
  
  unique_genes <- unique(use_markers_df$features.plot)
  
  # 2. 提取表达量数据
  # 我们调用Seurat内置的DotPlot函数，但不直接画图，而是获取其底层的计算数据
  p_base <- Seurat::DotPlot(obj, features = unique_genes, group.by = group.by)
  use_data <- p_base$data
  
  # 3. 将表达数据与我们自定义的分组元数据进行合并
  use_data <- merge(use_data, use_markers_df, by = "features.plot")
  
  # 4. 设置因子水平 (factor levels) 以严格保持用户输入的顺序
  use_data$features.plot <- factor(use_data$features.plot, levels = unique_genes)
  use_data$type <- factor(use_data$type, levels = group_names)
  
  # 5. 使用 ggplot2 重新构建分面图表
  p_final <- ggplot2::ggplot(use_data, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
    scale_size(range = size_range) + 
    scale_color_gradient2(low = colors[1], mid = colors[2], high = colors[3]) + 
    facet_wrap(~ type, scales = "free_x", nrow = 1) +  # 根据type分面，允许X轴自由缩放
    theme_classic(base_size = base_text_size) +
    theme(
      text = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      legend.title = element_text(size = base_text_size - 2), 
      strip.text = element_text(size = base_text_size - 1, face = "bold"), # 分面标签加粗
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1), # 分面框样式
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # 给每个分面添加边框
    ) +
    labs(x = "", y = "", size = "Percentage\nExpressed", color = "Scaled\nExpression")
  
  # 6. 可选：保存图表
  if (!is.null(save_prefix)) {
    ggsave(paste0(save_prefix, ".png"), plot = p_final, width = width, height = height)
    ggsave(paste0(save_prefix, ".pdf"), plot = p_final, width = width, height = height)
    message("图表已成功保存，文件前缀为: ", save_prefix)
  }
  
  return(p_final)
}



# # Seurat 按组分面气泡图 (Grouped DotPlot)

# 这是一个用于生成清晰、分面的单细胞RNA-seq气泡图（DotPlot）的R函数。与Seurat默认绘制的一长串连续基因不同，该函数允许您**按生物学功能、通路或特征基因集对基因进行分组**，使最终的可视化结果结构更清晰，更易于解读。

# ## 依赖包
# * `Seurat`
# * `ggplot2`

# ## 使用说明与示例

# 请先在您的 R 环境中加载该函数，然后将您感兴趣的基因整理为一个**命名列表 (Named List)**。列表的名称将作为图表顶部的分面标题。

# ```R
# # 1. 加载函数 (假设文件名为 GroupedDotPlot.R)
# source("GroupedDotPlot.R")

# # 2. 定义您的基因分组（使用命名列表）
# # 列表的名称 (如 "VM_structural_genes") 会自动成为分面标题
# my_gene_signatures <- list(
#   "VM_structural_genes" = c("CDH5", "KDR", "PECAM1", "CLDN5", "TEK", "CD34"),
#   "VM_regulatory_genes" = c("PDGFB", "ANGPT2", "VEGFA", "MMP9", "FN1", "COL1A1", "LOX"),
#   "Immunosuppressive_genes" = c("IL10", "TGFB1", "ARG1", "CD163", "MRC1", "CD274")
# )

# # 3. 生成并绘制图表
# # 注意：'seurat_object' 需要替换为您自己环境中的 Seurat 对象名称
# p <- GroupedDotPlot(
#   obj = seurat_object,
#   gene_groups = my_gene_signatures,
#   group.by = "celltype_myeloid",             # 映射到Y轴的细胞类型元数据列
#   colors = c("lightgrey", "#ffc7c7", "red"), # 自定义低中高表达量的颜色
#   save_prefix = "Dotplot_VM_genes_myeloid",  # 自动将结果保存为同名的 .png 和 .pdf 文件
#   width = 9,
#   height = 8
# )

# # 4. 在控制台/画板中显示图片
# print(p)





