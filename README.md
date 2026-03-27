# R Scripts: Visualization Results

本仓库收录了一些常用的 R 绘图函数及其输出示例，便于快速预览效果。

---

## 📊 Sankey_bar

Sankey 图，展示流向和分布。
```r
options(repr.plot.width=8, repr.plot.height=6)

# 绘制stage在细胞亚型中的流动分布
Sankey_bar(obj, 
             x_var = "sample",
             y_var = "stage",
             color_palette = cor_stage,
             output_prefix = "percent_stage_sub", width = 8, height = 6)
```


<p align="center">
  <img src="https://github.com/user-attachments/assets/2fedd35a-d8c4-4f6f-a549-490fe17bc7eb" width="500" />
</p>

---

## 📊 Sankey_bar_split

Sankey—split 图，展示分组流向和分布。

```r
Sankey_bar_split(obj,
         x_var = "maincelltype_0924",
         y_var = "day",
         group_var = "monkey",
         group_filter = c("M1","M2","M3","M4","M5"),
         output_prefix = "monkey_split",
         combine_plot = TRUE,
         show_labels = TRUE,
         label_threshold = 3)
```

         
<p align="center">
  <img width="668" height="841" alt="image" src="https://github.com/user-attachments/assets/42eb23e8-9bd4-489f-ab0c-1d06355eaef5" />

</p>


## 🌋 plot_volcano

火山图，常用于差异分析结果可视化。
```r
plot_volcano(
  data,
  output_file = "LUSC_tumor_volcano2.pdf",
  logFC_cutoff = 0.5,
  pval_cutoff = 0.05,
  highlight_genes = 'CLEC2B',
  label_top = 5,
  xlim = NULL,
  ylim = c(0, 350),
  gene_col = "X",  # 新增参数: 指定基因列
  width = 6, height = 6
)
```
<p align="center">
  <img src="https://github.com/user-attachments/assets/d18e3baf-518c-4580-b606-9297cb73c472" width="500" />
</p>

---

## 🧩 stacked_barplot

堆叠柱状图，展示类别组成比例。

```r
options(repr.plot.width=10, repr.plot.height=12)

stacked_barplot(
  obj = obj,
  group_x = "celltype_0818",
  group_y = "Time_category",
  file_name = "celltype_0818_Time_category",
  custom_colors = time_Colors,
  width = 10,
  height = 12
)
```
<p align="center">
  <img src="https://github.com/user-attachments/assets/40e327b5-3cc5-4a7c-b291-d8b1d7ce8d84" width="500" />
</p>





---


## Seurat 按组分面气泡图 (Grouped DotPlot)

绘制清晰、分面的单细胞RNA-seq气泡图（DotPlot）的R函数。与Seurat默认绘制的一长串连续基因不同，该函数允许您**按生物学功能、通路或特征基因集对基因进行分组**，使最终的可视化结果结构更清晰，更易于解读。

## 依赖包
* `Seurat`
* `ggplot2`

#####  使用说明与示例

请先在您的 R 环境中加载该函数，然后将您感兴趣的基因整理为一个**命名列表 (Named List)**。列表的名称将作为图表顶部的分面标题。

```R
# 1. 加载函数 (假设文件名为 GroupedDotPlot.R)
source("GroupedDotPlot.R")

# 2. 定义您的基因分组（使用命名列表）
# 列表的名称 (如 "VM_structural_genes") 会自动成为分面标题
my_gene_signatures <- list(
  "VM_structural_genes" = c("CDH5", "KDR", "PECAM1", "CLDN5", "TEK", "CD34"),
  "VM_regulatory_genes" = c("PDGFB", "ANGPT2", "VEGFA", "MMP9", "FN1", "COL1A1", "LOX"),
  "Immunosuppressive_genes" = c("IL10", "TGFB1", "ARG1", "CD163", "MRC1", "CD274")
)

# 3. 生成并绘制图表
# 注意：'seurat_object' 需要替换为您自己环境中的 Seurat 对象名称
p <- GroupedDotPlot(
  obj = seurat_object,
  gene_groups = my_gene_signatures,
  group.by = "celltype_myeloid",             # 映射到Y轴的细胞类型元数据列
  colors = c("lightgrey", "#ffc7c7", "red"), # 自定义低中高表达量的颜色
  save_prefix = "Dotplot_VM_genes_myeloid",  # 自动将结果保存为同名的 .png 和 .pdf 文件
  width = 9,
  height = 8
)

# 4. 在控制台/画板中显示图片
print(p)

<p align="center">
  <img width="2700" height="2100" alt="Dotplot_VM_genes_myeloid" src="https://github.com/user-attachments/assets/b604802d-2276-41de-990e-ad42e75ef541" />
</p>


---

## 📈 cluster_proportion_trend

按时间或条件的 cluster 比例变化趋势。

<p align="center">
  <img src="https://github.com/user-attachments/assets/982f0957-1ac8-496f-81bd-d19228c22485" width="400" />
</p>

---

## 📈 cluster_day_percent_by_monkey

统计每个猴子（monkey）在各时间点（day）中不同细胞群（group）的比例并绘图

```r
p <- cluster_day_percent_by_monkey_R(
  obj,
  group_col = "kc_1020",
  day_col = "day",
  monkey_col = "monkey",
  day_levels = c("1","3","7","14","21"),
  prefix = "kc1020_allMonkeys"
)
```

<p align="center">
  <img width="1648" height="985" alt="image" src="https://github.com/user-attachments/assets/a5f524ec-b55e-4c99-8305-d56e9ad7fe0f" />
</p>



