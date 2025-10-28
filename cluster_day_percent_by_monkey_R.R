# ============================================================
# Title: cluster_day_percent_by_monkey_R
# Author: Bolin & ChatGPT
# Date: 2025-10-27
# ============================================================

#' @title 统计每个猴子在各时间点（day）中不同细胞群（group）的比例并绘图
#'
#' @description
#' 本函数用于从 Seurat 对象中提取指定分组列（例如细胞亚群）、时间列（day）和猴子编号列，
#' 计算每个猴子在每个时间点下，各 group 所占的百分比（基于该猴子在该时间点的总细胞数），
#' 并生成按 group facet 的折线图。
#'
#' @param seurat_obj Seurat 对象，必须包含 @meta.data。
#' @param group_col 字符串，meta.data 中表示细胞类型或亚群的列名（例如 `"kc_1020"`）。
#' @param day_col 字符串，meta.data 中表示时间点的列名（例如 `"day"`）。
#' @param monkey_col 字符串，meta.data 中表示猴子编号的列名（例如 `"monkey"`）。
#' @param monkey_select 可选字符向量，仅分析指定猴子（如 `c("M1","M4")`），默认 `NULL` 表示全部。
#' @param day_levels 字符向量，设定时间点顺序（例如 `c("1","3","7","14","21")`）。
#' @param prefix 输出文件名前缀，例如 `"kc1020_by_monkey"`。
#'
#' @return 返回一个列表：
#' \itemize{
#'   \item plot: ggplot 对象，可继续自定义。
#'   \item data: 含比例信息的数据框 summary_full。
#' }
#'
#' @details
#' 计算逻辑如下：
#' 1. 提取指定列（group, day, monkey）；
#' 2. 统计每个 (day, group, monkey) 的细胞数；
#' 3. 在每个 (day, monkey) 内，将各 group 的计数除以该猴子该天的总细胞数，得出百分比；
#' 4. 自动补齐缺失组合（count=0），保证绘图面板完整；
#' 5. 自动识别 group_col 的 levels 顺序（若为因子）。
#'
#' @examples
#' \dontrun{
#' res <- cluster_day_percent_by_monkey_R(
#'   obj,
#'   group_col = "kc_1020",
#'   day_col = "day",
#'   monkey_col = "monkey",
#'   day_levels = c("1","3","7","14","21"),
#'   prefix = "kc1020_allMonkeys"
#' )
#' }
#'
#' @export
# ============================================================

cluster_day_percent_by_monkey_R <- function(
  seurat_obj,
  group_col = "kc_1020",
  day_col = "day",
  monkey_col = "monkey",
  monkey_select = NULL,
  day_levels = c("1", "3", "7", "14", "21"),
  prefix = "kc1020_by_monkey"
) {
  # ---- 输入检查 ----
  if (is.null(seurat_obj)) stop("❌ seurat_obj 为 NULL")
  meta_ok <- FALSE
  try({
    md <- seurat_obj@meta.data
    meta_ok <- !is.null(md)
  }, silent = TRUE)
  if (!meta_ok) stop("❌ seurat_obj 不包含可用的 @meta.data (Seurat 对象或 slot 不存在)。")

  meta <- seurat_obj@meta.data
  if (!(group_col %in% colnames(meta))) stop(paste0("❌ '", group_col, "' 不在 meta.data 中"))
  if (!(day_col %in% colnames(meta))) stop(paste0("❌ '", day_col, "' 不在 meta.data 中"))
  if (!(monkey_col %in% colnames(meta))) stop(paste0("❌ '", monkey_col, "' 不在 meta.data 中"))

  # ---- 提取并规范化数据 ----
  df <- meta[, c(group_col, day_col, monkey_col), drop = FALSE]
  colnames(df) <- c("group", "day", "monkey")
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # 自动继承 group 因子顺序
  if (is.factor(seurat_obj[[group_col]][, 1])) {
    group_levels <- levels(seurat_obj[[group_col]][, 1])
    df$group <- factor(df$group, levels = group_levels, ordered = TRUE)
  } else {
    df$group <- as.character(df$group)
  }

  df$day    <- factor(df$day, levels = as.character(day_levels), ordered = TRUE)
  df$monkey <- as.character(df$monkey)

  # ---- 统计 (day, group, monkey) 计数 ----
  summary_df <- df |>
    dplyr::group_by(day, group, monkey) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")

  # 过滤 monkey_select
  if (!is.null(monkey_select)) {
    monkey_select <- as.character(monkey_select)
    summary_df <- summary_df |>
      dplyr::filter(monkey %in% monkey_select)
    if (nrow(summary_df) == 0) stop("❌ 筛选后的数据为空，请检查 monkey_select 的取值。")
  }

  # ---- 计算每个 (day, monkey) 的总数并求比例 ----
  summary_df <- summary_df |>
    dplyr::group_by(day, monkey) |>
    dplyr::mutate(total_in_day_monkey = sum(count)) |>
    dplyr::ungroup() |>
    dplyr::mutate(proportion = ifelse(total_in_day_monkey > 0,
                                      count / total_in_day_monkey * 100, 0))

  # ---- 补齐缺失组合 ----
  all_days   <- levels(df$day)
  all_groups <- if (exists("group_levels")) group_levels else sort(unique(df$group))
  all_monkey <- if (is.null(monkey_select)) sort(unique(df$monkey)) else sort(unique(monkey_select))

  full_grid <- expand.grid(day = all_days, group = all_groups, monkey = all_monkey, stringsAsFactors = FALSE)
  full_grid$day <- factor(full_grid$day, levels = all_days, ordered = TRUE)

  summary_full <- dplyr::left_join(full_grid, summary_df, by = c("day", "group", "monkey")) |>
    dplyr::mutate(
      count = ifelse(is.na(count), 0, count),
      total_in_day_monkey = ifelse(is.na(total_in_day_monkey), 0, total_in_day_monkey)
    ) |>
    dplyr::group_by(day, monkey) |>
    dplyr::mutate(total_in_day_monkey = ifelse(total_in_day_monkey == 0 & sum(count) > 0,
                                               sum(count), total_in_day_monkey)) |>
    dplyr::ungroup() |>
    dplyr::mutate(proportion = ifelse(total_in_day_monkey > 0,
                                      count / total_in_day_monkey * 100, 0))

  if (exists("group_levels")) {
    summary_full$group <- factor(summary_full$group, levels = group_levels, ordered = TRUE)
  }

  message("✅ 数据统计完成，预览前几行：")
  print(head(summary_full))
  utils::write.csv(summary_full, paste0(prefix, "_summary_full.csv"), row.names = FALSE)

  # ---- 绘图 ----
  p <- ggplot2::ggplot(summary_full, ggplot2::aes(x = day, y = proportion, group = monkey, color = monkey)) +
    ggplot2::geom_line(ggplot2::aes(linetype = monkey), size = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~group, scales = "fixed", ncol = 3) +
    ggplot2::labs(
      x = "Day", y = "Proportion (%)",
      title = paste0("Cluster proportions over days by monkey (", group_col, ")")
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(size = 11),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  out_file <- paste0(prefix, "_", group_col, "_by_monkey_percent.png")
  ggplot2::ggsave(out_file, p, width = 20, height = 12, dpi = 300)
  message("✅ 图片已保存为: ", out_file)

  invisible(list(plot = p, data = summary_full))
}


# source("cluster_day_percent_by_monkey_R.R")

# res_all <- cluster_day_percent_by_monkey_R(
#   obj,
#   group_col = "kc_1020",
#   day_col = "day",
#   monkey_col = "monkey",
#   day_levels = c("1","3","7","14","21"),
#   prefix = "kc1020_allMonkeys"
# )
