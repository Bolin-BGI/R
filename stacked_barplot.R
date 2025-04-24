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
stacked_barplot <- function(obj, group_x, group_y, file_name='stacked_barplot.png', width = 12, height = 8, dpi = 300) {
  
  # 提取元数据
  meta <- obj@meta.data %>%
    select(!!sym(group_x), !!sym(group_y)) %>%
    filter(!is.na(!!sym(group_x)))  # 移除缺失值
  
  # 统计细胞数量
  counts <- meta %>%
    group_by(!!sym(group_x), !!sym(group_y)) %>%
    tally(name = "count") %>%
    ungroup()
  
  # 转换为比例
  prop <- counts %>%
    group_by(!!sym(group_x)) %>%
    mutate(percent = (count / sum(count)) * 100)
  
  # 定义颜色方案
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
  
  # 确保颜色数量与分组匹配
  unique_groups <- unique(prop[[group_y]])
  colors <- colorlist[seq_along(unique_groups)]
  
  # 绘图
  p <- ggplot(prop, 
              aes(x = !!sym(group_x), 
                  y = percent, 
                  fill = !!sym(group_y))) +
    geom_col() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),          # 移除网格线
      axis.text.x = element_text(
        angle = 45, 
        hjust = 1, 
        face = "bold"                       # 横坐标标签加粗
      ),
      axis.title = element_text(face = "bold"),  # 坐标轴标题加粗
      plot.title = element_text(face = "bold")   # 图表标题加粗（若有）
    ) +
    labs(
      x = "Cell Type", 
      y = "Percentage", 
      fill = "Group"
    ) +
    scale_fill_manual(values = colors)  # 应用颜色方案
  
  # 添加标签
  p <- p + geom_text(
    aes(label = sprintf("%.1f%%", percent)),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "white",
    fontface = "bold"  # 标签文字加粗
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
  
  # 返回 ggplot 对象以便后续调整
  return(p)
}