# R Scripts: Visualization Results

本仓库收录了一些常用的 R 绘图函数及其输出示例，便于快速预览效果。

---

## 📊 Sankey_bar

Sankey 图，展示流向和分布。

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

<p align="center">
  <img src="https://github.com/user-attachments/assets/69fcccbc-92a4-4a83-b9a6-5b62e123bc14" width="300" />
  <img src="https://github.com/user-attachments/assets/d18e3baf-518c-4580-b606-9297cb73c472" width="500" />
</p>

---

## 🧩 stacked_barplot

堆叠柱状图，展示类别组成比例。

<p align="center">
  <img src="https://github.com/user-attachments/assets/40e327b5-3cc5-4a7c-b291-d8b1d7ce8d84" width="500" />
</p>

---

## 📈 cluster_proportion_trend

按时间或条件的 cluster 比例变化趋势。

<p align="center">
  <img src="https://github.com/user-attachments/assets/982f0957-1ac8-496f-81bd-d19228c22485" width="400" />
</p>
