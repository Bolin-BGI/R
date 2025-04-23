plot_volcano <- function(
  data, 
  output_file = "volcano.pdf",
  logFC_cutoff = 0.5, 
  pval_cutoff = 0.05,
  highlight_genes = NULL,
  label_top = 5,
  xlim = NULL,
  ylim = c(0, 350),
  # 新增参数：指定基因列名
  gene_col = "gene"
) {
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  
  # 数据校验
  if (!gene_col %in% colnames(data)) {
    stop(paste0("基因列 '", gene_col, "' 不存在于数据中，请检查参数gene_col设置"))
  }
  
  # 数据预处理
  df <- data %>%
    rename(
      logFC = avg_log2FC,
      FDR = p_val_adj
    ) %>%
    mutate(
      # 处理极端p值
      FDR = ifelse(FDR <= 0, 1e-323, FDR),
      # 创建显著性标签
      significance = case_when(
        FDR < pval_cutoff & logFC > logFC_cutoff ~ "Up",
        FDR < pval_cutoff & logFC < -logFC_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # 自动标签逻辑
  auto_labels <- df %>%
    filter(significance != "NS") %>%
    group_by(significance) %>%
    arrange(FDR, desc(abs(logFC))) %>%
    slice_head(n = label_top)
  
  # 手动标签逻辑
  if (!is.null(highlight_genes)) {
    manual_labels <- df %>% 
      filter(!!sym(gene_col) %in% highlight_genes)
    label_data <- bind_rows(auto_labels, manual_labels) %>%
      distinct(!!sym(gene_col), .keep_all = TRUE)
  } else {
    label_data <- auto_labels
  }
  
  # 颜色映射
  color_palette <- c(Up = "#E41A1C", Down = "#377EB8", NS = "grey60")
  
  # 创建基础绘图对象
  p <- ggplot(df, aes(logFC, -log10(FDR))) +
    geom_point(aes(color = significance), 
               size = 2, alpha = 0.7) +
    scale_color_manual(values = color_palette) +
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), 
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(pval_cutoff), 
               linetype = "dashed", color = "grey40") +
    labs(x = expression(log[2]("Fold Change")),
         y = expression(-log[10]("FDR")),
         color = "Significance") +
    theme_classic(base_size = 12) +
    theme(legend.position = "top")
  
  # 添加标签
  if (nrow(label_data) > 0) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = label_data,
        aes_string(label = gene_col),  # 改用aes_string
        size = 3,
        box.padding = 0.3,
        segment.color = "grey50",
        max.overlaps = Inf
      )
  }
  
  # 坐标轴调整
  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(ylim)) p <- p + ylim(ylim)
  
  # 保存图形
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
}