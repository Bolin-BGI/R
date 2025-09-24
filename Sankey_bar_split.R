# =========================================================
# 函数名称: Sankey_bar_split
# 功能: 绘制分组变量的堆叠流图 (Sankey-like barplot)
# 适用对象: Seurat 对象的 meta.data 或普通 data.frame
#
# 参数说明:
#   obj            : Seurat 对象 或 data.frame (必须包含 x_var, y_var, group_var 三列)
#   x_var          : 字符串，表示横坐标分组变量名（常用于细胞类型）
#   y_var          : 字符串，表示纵坐标分组变量名（常用于时间点）
#   group_var      : 字符串，表示分面变量名（常用于样本、个体、批次）
#   group_filter   : 可选，向量，指定要保留的 group_var 取值；NULL 表示不过滤
#   color_palette  : 可选，自定义颜色向量，顺序需与 x_var 的 levels 对应
#   output_prefix  : 输出文件名前缀
#   width, height  : 输出图像的宽高
#   combine_plot   : TRUE 表示将所有分组画在同一张图（facet_wrap）
#   use_percent    : TRUE 计算百分比，FALSE 显示原始计数
#   exclude_zero_y : TRUE 过滤掉 y_var == 0 的行，适合去掉 Day0
#   show_labels    : TRUE 在柱子中显示百分比或数值
#   label_threshold: 标签显示阈值（低于该百分比不显示）
#
# 使用示例:
#   Sankey_bar_split(obj,
#              x_var = "maincelltype_0924",
#              y_var = "day",
#              group_var = "monkey",
#              group_filter = c("M1","M2","M3","M4","M5"),
#              output_prefix = "stage_composition",
#              combine_plot = TRUE,
#              show_labels = TRUE,
#              label_threshold = 3)
#
# =========================================================

Sankey_bar_split <- function(obj,
                       x_var = "celltype",
                       y_var = "time",
                       group_var = "sample",
                       group_filter = NULL,
                       color_palette = NULL,
                       output_prefix = "sankey_plot",
                       width = 14,
                       height = 8,
                       combine_plot = TRUE,
                       use_percent = TRUE,
                       exclude_zero_y = TRUE,
                       show_labels = TRUE,
                       label_threshold = 3) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  
  # -------- 数据准备 --------
  # 如果是 Seurat 对象，提取 meta.data，否则直接用 data.frame
  if ("Seurat" %in% class(obj)) {
    df <- obj@meta.data %>% tibble::rownames_to_column("cell_id")
  } else {
    df <- obj %>% tibble::as_tibble()
    if (!"cell_id" %in% colnames(df)) df <- df %>% tibble::rowid_to_column("cell_id")
  }
  
  # 按 group_filter 筛选
  if (!is.null(group_filter)) {
    df <- df %>% dplyr::filter(.data[[group_var]] %in% group_filter)
  }
  
  # 可选：过滤掉 y_var == 0
  if (exclude_zero_y) {
    df <- df %>% dplyr::filter(.data[[y_var]] != 0)
  }
  
  # 统计每个组合的细胞数量
  df_count <- df %>%
    group_by(.data[[group_var]], .data[[y_var]], .data[[x_var]]) %>%
    summarise(count = n(), .groups = "drop")
  
  # 百分比或原始数值
  if (use_percent) {
    df_count <- df_count %>%
      group_by(.data[[group_var]], .data[[y_var]]) %>%
      mutate(value = count / sum(count) * 100) %>%
      ungroup()
  } else {
    df_count <- df_count %>% mutate(value = count)
  }
  
  # 保持 factor 顺序一致
  df_count[[x_var]] <- factor(df_count[[x_var]], levels = unique(df_count[[x_var]]))
  df_count[[y_var]] <- factor(df_count[[y_var]], levels = sort(unique(df_count[[y_var]])))
  
  # -------- 绘图子函数 --------
  build_plot <- function(plot_df) {
    
    # 提前生成标签
    plot_df <- plot_df %>% 
      mutate(lab = if (use_percent) sprintf("%.1f%%", value) else as.character(value))
    
    p <- ggplot(plot_df, aes(x = .data[[y_var]], y = value, fill = .data[[x_var]])) +
      # 用透明柱确定堆叠高度，供 label 定位
      geom_col(position = position_stack(), width = 0.5, alpha = 0) +
      # 绘制 alluvial 流
      geom_flow(aes(stratum = .data[[x_var]], alluvium = .data[[x_var]]),
                width = 0.5, alpha = 0.7, curve_type = "cubic") +
      geom_stratum(aes(stratum = .data[[x_var]]),
                   width = 0.5, color = "black", size = 0.2) +
      scale_fill_manual(values = if (is.null(color_palette)) {
        pal0 <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941",
                  "#006FA6","#A30059","#FFE4E1","#0000A6","#63FFAC")
        rep(pal0, length.out = length(unique(plot_df[[x_var]])))
      } else color_palette, drop = FALSE) +
      labs(x = y_var, y = ifelse(use_percent, "Percentage (%)", "Count"), fill = x_var) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 12))
    
    # 添加标签
    if (show_labels) {
      label_data <- plot_df %>% filter(value > label_threshold)
      if (nrow(label_data) > 0) {
        p <- p + geom_text(data = label_data,
                           aes(x = .data[[y_var]], y = value, label = lab, group = .data[[x_var]]),
                           position = position_stack(vjust = 0.5),
                           size = 2, inherit.aes = FALSE)
      }
    }
    return(p)
  }
  
  # -------- 输出部分 --------
  if (combine_plot) {
    # 合并绘制 (facet_wrap)
    p <- build_plot(df_count) + facet_wrap(as.formula(paste("~", group_var)))
    ggsave(paste0(output_prefix, "_combined.png"), p,
           width = width, height = height, bg = "white", dpi = 300)
    print(p)
  } else {
    # 按分组单独输出
    for (g in unique(df_count[[group_var]])) {
      df_g <- df_count %>% filter(.data[[group_var]] == g)
      p <- build_plot(df_g) + ggtitle(paste(group_var, g))
      ggsave(paste0(output_prefix, "_", g, ".png"), p,
             width = width, height = height, bg = "white", dpi = 300)
      print(p)
    }
  }
}
