# 定义 Sankey_bar 函数
Sankey_bar <- function(obj, 
                       x_var = "sub",  # 横坐标变量名（默认为 "sub"）
                       y_var = "stage",  # 纵坐标变量名（默认为 "stage"）
                       color_palette = NULL,  # 颜色调色板（可选）
                       output_prefix = "stage_composition",  # 输出文件名前缀
                       width = 8,  # 图片宽度（单位：英寸）
                       height = 10) {  # 图片高度（单位：英寸）
  
  # 加载必要的库
  library(ggalluvial)  # 用于绘制 Sankey 图
  library(data.table)  # 数据处理
  
  # 提供的颜色列表（预定义的大量颜色）
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
  
  # 数据过滤与格式转换
  # 筛选出 x_var 和 y_var 列均不为空的数据
  filtered_data <- obj@meta.data %>% dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) 
  
  # 核心比例计算逻辑
  # 计算每个 x_var 类别中 y_var 的占比
  percentage <- as.data.frame(table(filtered_data[[x_var]], filtered_data[[y_var]])) %>%
    setNames(c(x_var, y_var, "count")) %>%
    group_by(.data[[x_var]]) %>%  # 按 x_var 分组
    mutate(total = sum(count),  # 计算每组总和
           percentage = count / total * 100)  # 计算百分比
  
  # 因子化处理，并确保 x_var 的顺序与 obj@meta.data 的因子顺序一致
  percentage[[x_var]] <- factor(percentage[[x_var]], 
                               levels = levels(obj@meta.data[[x_var]]))  # 保持 x_var 的因子顺序
  percentage[[y_var]] <- factor(percentage[[y_var]], 
                               levels = unique(filtered_data[[y_var]]))  # 确保 y_var 的因子顺序
  
  # 如果 color_palette 未指定，则从 colorlist 中按顺序取颜色
  if (is.null(color_palette)) {
    num_colors <- length(unique(filtered_data[[y_var]]))  # 获取 y_var 的唯一值数量
    color_indices <- rep(1:length(colorlist), length.out = num_colors)  # 循环使用 colorlist
    color_palette <- colorlist[color_indices]  # 从 colorlist 中按顺序取颜色
    
    # 将颜色向量命名为 y_var 的唯一值
    names(color_palette) <- unique(filtered_data[[y_var]])
  }
  
  # 可视化层（保持流式图核心元素）
  p <- ggplot(percentage, 
              aes(x = .data[[x_var]], y = percentage,  # x 轴为 x_var，y 轴为百分比
                  stratum = .data[[y_var]],  # 层级变量为 y_var
                  alluvium = .data[[y_var]],  # 流线变量为 y_var
                  fill = .data[[y_var]])) +  # 填充颜色由 y_var 决定
    geom_flow(width = 0.5, curve_type = "linear", alpha = 0.5) +  # 绘制流线
    geom_stratum(width = 0.5, alpha = 1, color = NA) +  # 绘制层级块
    geom_alluvium(width = 0.5, curve_type = "linear",  # 绘制流线
                  color = "#f4e2de", fill = NA, size = 0.8) +
    scale_fill_manual(values = color_palette) +  # 设置填充颜色
    labs(x = "Celltype", y = "Percentage (%)", fill = y_var) +  # 设置轴标签和图例标题
    theme(
      panel.grid = element_blank(),  # 移除网格线
      panel.background = element_blank(),  # 移除背景填充
      plot.background = element_rect(fill = "white"),  # 设置绘图背景为白色
      panel.border = element_rect(fill = NA, color = "black", size = 1),  # 添加边框
      axis.text.x = element_text(angle = 0, hjust = 1, color = "black", size = 13),  # 设置 x 轴文本样式
      axis.text.y = element_text(color = "black", size = 15),  # 设置 y 轴文本样式
      axis.title = element_text(size = 16),  # 设置轴标题字体大小
      legend.text = element_text(size = 10)  # 设置图例文字大小
    )
  
  # 输出文件
  ggsave(paste0(output_prefix, ".png"), p, width = width, height = height, bg='white', dpi = 300)  # 保存为 PNG 文件
  ggsave(paste0(output_prefix, ".pdf"), p, width = width, height = height)  # 保存为 PDF 文件
  
  return(p)  # 返回生成的 ggplot 对象
}
