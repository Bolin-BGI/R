#' stacked_barplot: 绘制细胞类型分组的堆积柱状图
#'
#' @param obj Seurat 对象，包含元数据（需包含细胞类型和分组列）
#' @param group_x 字符串，指定元数据中细胞类型列的名称（如 "kc_celltype"）
#' @param group_y 字符串，指定元数据中分组列的名称（如 "monkey"）
#' @param file_name 字符串，保存的图片文件名（如 "stacked_plot.png"）
#' @param width 数值，图片宽度（英寸）
#' @param height 数值，图片高度（英寸）
#' @param dpi 数值，分辨率（默认 300）
#'
#' @return 返回 ggplot 对象（可进一步调整或查看）
#'
#' @examples
#' stacked_barplot(
#'   obj = your_seurat_object,
#'   group_x = "kc_celltype",
#'   group_y = "monkey",
#'   file_name = "output.png",
#'   width = 12,
#'   height = 8
#' )

Sankey_bar <- function(obj, 
                       x_var = "sub", 
                       y_var = "stage", 
                       color_palette = NULL, 
                       output_prefix = "stage_composition", 
                       width = 8, 
                       height = 10) {
  
  library(ggalluvial)
  library(data.table)
  library(dplyr)
  library(tidyr)
  
  colorlist <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", 
                 "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
                 "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
                 "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
                 "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
                 "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
                 "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
                 "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
                 "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
                 "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
                 "#191970", "#000080", "#6495ED", "#1E90FF", "#00BFFF",
                 "#00FFFF", "#FF1493", "#FF00FF", "#A020F0", "#63B8FF",
                 "#008B8B", "#54FF9F", "#00FF00", "#76EE00", "#FFF68F")
  
  # 数据过滤
  filtered_data <- obj@meta.data %>% 
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]]), 
                  .data[[x_var]] != "", .data[[y_var]] != "")
  
  # 调试：检查过滤后数据
  print("Filtered data dimensions:")
  print(dim(filtered_data))
  # print("Head of filtered data:")
  # print(head(filtered_data))
  
  # 计算比例
  percentage <- as.data.frame(table(filtered_data[[x_var]], filtered_data[[y_var]])) %>%
    setNames(c(x_var, y_var, "count")) %>%
    complete(.data[[x_var]], .data[[y_var]], fill = list(count = 0)) %>%
    group_by(.data[[x_var]]) %>%
    mutate(total = sum(count), percentage = count / total * 100)
  
  # 调试：检查 percentage 数据框
  # print("Percentage data frame:")
  # print(head(percentage))
  # print(str(percentage))
  
  # 因子化
  percentage[[x_var]] <- factor(percentage[[x_var]], 
                               levels = unique(percentage[[x_var]]))
  percentage[[y_var]] <- factor(percentage[[y_var]], 
                               levels = unique(percentage[[y_var]]))
  
  # 设置颜色
  if (is.null(color_palette)) {
    num_colors <- length(unique(filtered_data[[y_var]]))
    color_indices <- rep(1:length(colorlist), length.out = num_colors)
    color_palette <- colorlist[color_indices]
    names(color_palette) <- unique(filtered_data[[y_var]])
  }
  
  # 绘图
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
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 18),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 13),
      plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")  # 调整边距: 顶部、右边、底部、左边
    )
  
  # 保存图片
  ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, bg = "white", dpi = 300)
  ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)
  
  return(p)
}






# stacked_barplot <- function(obj, group_x, group_y, file_name='stacked_barplot.png', width = 12, height = 8, dpi = 300) {
  
#   # 提取元数据
#   meta <- obj@meta.data %>%
#     select(!!sym(group_x), !!sym(group_y)) %>%
#     filter(!is.na(!!sym(group_x)))  # 移除缺失值
  
#   # 统计细胞数量
#   counts <- meta %>%
#     group_by(!!sym(group_x), !!sym(group_y)) %>%
#     tally(name = "count") %>%
#     ungroup()
  
#   # 转换为比例
#   prop <- counts %>%
#     group_by(!!sym(group_x)) %>%
#     mutate(percent = (count / sum(count)) * 100)
  
#   # 定义颜色方案
#   colorlist <- c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", 
#                  "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
#                  "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
#                  "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
#                  "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
#                  "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
#                  "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
#                  "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
#                  "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
#                  "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
#                  "#191970", "#000080", "#6495ED", "#1E90FF", "#00BFFF",
#                  "#00FFFF", "#FF1493", "#FF00FF", "#A020F0", "#63B8FF",
#                  "#008B8B", "#54FF9F", "#00FF00", "#76EE00", "#FFF68F")
  
#   # 确保颜色数量与分组匹配
#   unique_groups <- unique(prop[[group_y]])
#   colors <- colorlist[seq_along(unique_groups)]
  
#   # 绘图
#   p <- ggplot(prop, 
#               aes(x = !!sym(group_x), 
#                   y = percent, 
#                   fill = !!sym(group_y))) +
#     geom_col() +
#     theme_minimal() +
#     theme(
#       panel.grid = element_blank(),          # 移除网格线
#       axis.text.x = element_text(
#         angle = 45, 
#         hjust = 1, 
#         face = "bold"                       # 横坐标标签加粗
#       ),
#       axis.title = element_text(face = "bold"),  # 坐标轴标题加粗
#       plot.title = element_text(face = "bold")   # 图表标题加粗（若有）
#     ) +
#     labs(
#       x = "Cell Type", 
#       y = "Percentage", 
#       fill = "Group"
#     ) +
#     scale_fill_manual(values = colors)  # 应用颜色方案
  
#   # 添加标签
#   p <- p + geom_text(
#     aes(label = sprintf("%.1f%%", percent)),
#     position = position_stack(vjust = 0.5),
#     size = 3,
#     color = "white",
#     fontface = "bold"  # 标签文字加粗
#   )
  
#   # 保存图片
#   ggsave(
#     filename = file_name,
#     plot = p,
#     width = width,
#     height = height,
#     dpi = dpi,
#     bg = "white"
#   )
  
#   # 返回 ggplot 对象以便后续调整
#   return(p)
# }
