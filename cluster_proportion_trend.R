##############################################################################
# 使用示例：
# 加载调色板包（首次使用需安装）
# install.packages(c("ggsci", "RColorBrewer", "viridis"))
# 
# cluster_proportion_trend(
#   seurat_obj = your_seurat_object,
#   cluster_col = "seurat_clusters",
#   group_col = "stage",
#   palette_name = "D3", # 尝试更换调色板："Set3","Viridis"
#   prefix = "my_analysis"
# )
##############################################################################


#' @title 绘制细胞群比例随分组变量变化趋势（改进调色板版）
#' @description 可视化指定聚类列在不同分组条件下的比例变化趋势，支持多种调色板
#' @param seurat_obj Seurat 对象
#' @param cluster_col 表示细胞群的列名 (e.g. "seurat_clusters")
#' @param group_col 表示分组变量的列名 (e.g. "stage")
#' @param palette_name 调色板名称，支持："D3","Set3","Paired","Spectral","Viridis" (default: "D3")
#' @param prefix 输出图片前缀 (default: "")
#' @param show_plot 是否显示绘图结果 (default: TRUE)
cluster_proportion_trend <- function(seurat_obj, 
                                    cluster_col, 
                                    group_col, 
                                    palette_name = "D3",
                                    prefix = "",
                                    show_plot = TRUE) {
  # 加载必要包
  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(scales)
  require(ggsci)    # 新增专业调色板
  require(RColorBrewer)
  require(viridis)
  
  # 提取元数据
  metadata <- seurat_obj@meta.data %>% 
                tibble::rownames_to_column("cell") %>% 
                select(all_of(c(cluster_col, group_col))) %>% 
                mutate(!!group_col := factor(.data[[group_col]]))  # 确保分组列为因子
  
  # 计算比例
  prop_data <- metadata %>%
                group_by(.data[[group_col]], .data[[cluster_col]]) %>% 
                summarise(count = n(), .groups = "drop_last") %>% 
                mutate(proportion = count / sum(count) * 100) %>% 
                ungroup()
  
  # 获取调色板（新增多种专业调色板）
  n_clusters <- n_distinct(prop_data[[cluster_col]])
  color_palette <- switch(
    palette_name,
    "D3"       = pal_d3("category20")(n_clusters),
    "Set3"     = brewer.pal(n_clusters, "Set3"),
    "Paired"   = brewer.pal(n_clusters, "Paired"),
    "Spectral" = brewer.pal(n_clusters, "Spectral"),
    "Viridis"  = viridis(n_clusters),
    pal_d3("category20")(n_clusters) # 默认使用D3
  )
  
  # 绘图参数设置
  common_theme <- theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.background = element_blank())
  
  # 分面折线图
  p_facet <- ggplot(prop_data, aes(x = .data[[group_col]], 
y = proportion,
group = 1,
 color = .data[[cluster_col]])) +
                geom_line(linewidth = 1) +
                geom_point(size = 3) +
                facet_wrap(as.formula(paste("~", cluster_col)), ncol = 4) +
                scale_color_manual(values = color_palette) +
                labs(title = "Cluster Proportion by Group",
                     x = group_col,
                     y = "Proportion (%)") +
                common_theme +
                theme(legend.position = "none")
  
  # 组合折线图
  p_combined <- ggplot(prop_data, aes(x = .data[[group_col]], 
                                    y = proportion,
                                    group = .data[[cluster_col]],
                                    color = .data[[cluster_col]])) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = color_palette) +
    labs(title = "Combined Cluster Proportions",
         x = group_col,
         y = "Proportion (%)",
         color = cluster_col) +
    common_theme +
    guides(color = guide_legend(ncol = 2)) # 分两列显示图例
  
  # 保存图片
  if (prefix != "") {
    ggsave(paste0(prefix, "_facet_proportion.png"), 
           p_facet, width = 15, height = 10, dpi = 300)
    ggsave(paste0(prefix, "_combined_proportion.png"), 
           p_combined, width = 12, height = 8, dpi = 300)
  }
  
  # 显示图形
  if (show_plot) {
    print(p_facet)
    print(p_combined)
  }
  
  # 返回计算结果
  invisible(prop_data)
}

