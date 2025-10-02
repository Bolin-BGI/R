#' stacked_barplot: 绘制细胞类型分组的堆积柱状图
#'
#' @param obj Seurat 对象，包含元数据（需包含细胞类型和分组列）
#' @param group_x 字符串，指定元数据中细胞类型列的名称（如 "kc_celltype"）
#' @param group_y 字符串，指定元数据中分组列的名称（如 "monkey"）
#' @param file_name 字符串，保存的图片文件名（如 "stacked_plot.png"）
#' @param width 数值，图片宽度（英寸）
#' @param height 数值，图片高度（英寸）
#' @param dpi 数值，分辨率（默认 300）
#' @param custom_colors 向量，可选。指定每个分组的颜色，命名应与 group_y 的水平一致。
#'
#' @return 返回 ggplot 对象（可进一步调整或查看）
#'


# example
# stacked_barplot(
#   obj = obj,
#   group_x = "celltype_0818",
#   group_y = "Grade",
#   file_name = "output.png",
#   custom_colors = c(M1 = "#E64B35FF", M2 = "#4DBBD5FF",****),
#   width = 12,
#   height = 8
# )


stacked_barplot <- function(obj, 
                            group_x, 
                            group_y, 
                            file_name = 'stacked_barplot.png', 
                            width = 12, 
                            height = 8, 
                            dpi = 300, 
                            custom_colors = NULL) {
  
  # 提取元数据
  meta <- obj@meta.data %>%
    dplyr::select(!!sym(group_x), !!sym(group_y)) %>%
    dplyr::filter(!is.na(!!sym(group_x)))  # 移除缺失值
  
  # 统计细胞数量
  counts <- meta %>%
    group_by(!!sym(group_x), !!sym(group_y)) %>%
    tally(name = "count") %>%
    ungroup()
  
  # 转换为比例
  prop <- counts %>%
    group_by(!!sym(group_x)) %>%
    mutate(percent = (count / sum(count)) * 100)
  
  # 预定义颜色方案
  colorlist <- c("#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
                 "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", 
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
  
  # 确定颜色方案
  unique_groups <- unique(prop[[group_y]])
  if (!is.null(custom_colors)) {
    # 如果用户提供了 custom_colors，检查并使用
    if (!all(unique_groups %in% names(custom_colors))) {
      stop("custom_colors 必须包含所有分组的命名颜色")
    }
    colors <- custom_colors[unique_groups]
  } else {
    # 否则使用默认方案
    colors <- colorlist[seq_along(unique_groups)]
    names(colors) <- unique_groups
  }
  
  # 绘图
  p <- ggplot(prop, 
              aes(x = !!sym(group_x), 
                  y = percent, 
                  fill = !!sym(group_y))) +
    geom_col() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),          
      axis.text.x = element_text(
        angle = 45, 
        hjust = 0.6, 
        size = 12,
        face = "bold"                       
      ),
      axis.title = element_text(face = "bold"),  
      plot.title = element_text(face = "bold"),   
      plot.margin = margin(t = 5, r = 5, b = 10, l = 20, unit = "mm") 
    ) +
    labs(
      x = "Cell Type", 
      y = "Percentage", 
      fill = "Group"
    ) +
    scale_fill_manual(values = colors)  
  
  
    # 添加标签（只显示 >=1% 的部分）
    p <- p + geom_text(
      data = subset(prop, percent >= 1),  # 筛选出百分比 >=1 的行
      aes(label = sprintf("%.1f%%", percent),
          x = !!sym(group_x), 
          y = percent,
          fill = !!sym(group_y)),
      position = position_stack(vjust = 0.5),
      size = 4,
      color = "white",
      fontface = "bold"
    )
  
  # 保存图片
  ggsave(
    filename = file_name,
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  return(p)
}
