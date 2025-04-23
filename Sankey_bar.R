#' 绘制细胞类型-阶段占比流式图
#' @param obj Seurat对象
#' @param x_var 横坐标变量名（细胞类型列）
#' @param y_var 纵坐标变量名（阶段分组列）
#' @param color_palette 配色向量（需命名）
#' @param output_prefix 输出文件名前缀
#' @param width 图片宽度（单位：英寸，默认10）
#' @param height 图片高度（单位：英寸，默认12）
Sankey_bar <- function(obj, 
                                  x_var = "sub",
                                  y_var = "stage",
                                  color_palette = cor_stage,
                                  output_prefix = "stage_composition",
                                  width = 8,
                                  height = 10) {
    library(ggalluvial)
    library(data.table)
    library(ggsci)

  # 数据过滤与格式转换
  filtered_data <- obj@meta.data %>% dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) 
  
  # 核心比例计算逻辑（保持原有结构）
  percentage <- as.data.frame(table(filtered_data[[x_var]], filtered_data[[y_var]])) %>%
    setNames(c(x_var, y_var, "count")) %>%
    group_by(.data[[x_var]]) %>%  # 按细胞类型分组
    mutate(total = sum(count), 
           percentage = count / total * 100)
  
  # 因子化处理
  percentage[[x_var]] <- factor(percentage[[x_var]], 
                               levels = unique(filtered_data[[x_var]]))
  percentage[[y_var]] <- factor(percentage[[y_var]], 
                               levels = names(color_palette))
  
  # 可视化层（保持流式图核心元素）
  p <- ggplot(percentage, 
              aes(x = .data[[x_var]], y = percentage,
                  stratum = .data[[y_var]], 
                  alluvium = .data[[y_var]],
                  fill = .data[[y_var]])) +
    geom_flow(width = 0.5, curve_type = "linear", alpha = 0.5) +
    geom_stratum(width = 0.5, alpha = 1, color = NA) +
    geom_alluvium(width = 0.5, curve_type = "linear", 
                  color = "#f4e2de", fill = NA, size = 0.8) +
    scale_fill_manual(values = color_palette) +
    labs(x = "Celltype", y = "Percentage (%)", fill = y_var) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white"),
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 13),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 10)
    )
  
  # 输出文件
  ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, bg='white', dpi = 300)
  ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
  
  return(p)
}