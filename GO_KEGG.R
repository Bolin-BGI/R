library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)


# 读取CSV文件
marker <- read.csv("./kc_celltype_findmarkers.csv",row.names = 1)
colnames(marker)[1] <- "geneID" # 将第一列命名为'gene'

unique(marker$cluster)
# [1] "KC-migratory"         "Activated KC_spinous"
# [3] "TAC"                  "KC_basal_prolif"     
# [5] "KC-HS"                "KC-spinous"          
# [7] "KC_Stress-Responsive" "KC-basal"            
# [9] "KC_ACTA2"             "KC-IRS"              
# [11] "SAC"                  "KC-ORS

name = 'KC_Stress-Responsive' # 'KC-migratory'  'Activated KC_spinous' 'KC_basal_prolif' 'KC_Stress-Responsive'

# 筛选上调表达的基因
top_up <- marker %>% filter(cluster == name) %>% filter(avg_log2FC > 0.5 & p_val < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
  filter(!grepl("^ENSSS", gene) & !grepl("^ENSMFA", gene) & !grepl("^(RPL|RPS)", gene)) # %>% head(100)
nrow(top_up)

# 转换基因符号为 ENTREZID
gene_ID <- bitr(top_up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene_ID)

# ----------------------------------- GO ----------------------------------

# 设置p值和q值过滤阈值
pvalueFilter = 0.05  
qvalueFilter = 0.05

# 根据q值过滤阈值选择颜色
colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

# 设置显示的条目数
showNum = 10


# main
GO <- enrichGO(gene = gene_ID$ENTREZID,OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",  # 记得和输入基因的格式保持一致
               ont = "ALL", # BP(Biological Process)  MF(Molecular Function)  CC(Cellular Component) 
               pvalueCutoff = 0.5, 
               qvalueCutoff = 0.5,
               readable = TRUE)
print(paste0('The number of GO_items is : ',nrow(GO)))

df = as.data.frame(GO)
df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
write.csv(df_ed, file = "positive_GO_up.csv", row.names = FALSE)

if (nrow(GO) < 10) {
  showNum = nrow(GO)
}

# 如果有富集结果，生成气泡图
if (nrow(GO) != 0) {
  # 生成气泡图
  pdf(file = paste0(name, "_GO.pdf"), width = 12, height = 15)
  p1 <- dotplot(GO,x = "GeneRatio", color = colorSel, size = "Count", showCategory = showNum,label_format=150,split="ONTOLOGY") + # 以 ONTOLOGY 类型分开  
    facet_grid(ONTOLOGY~., scales = 'free') # 以 ONTOLOGY 类型分屏绘图
  print(p1)
  #dev.off()
  
  
  # 生成条形图
  #pdf(file = "GO_bar.pdf", width = 10, height = 8)
  p2 <- barplot(GO, x = "Count", color = colorSel, showCategory = showNum, label_format=150,split="ONTOLOGY") + 
    facet_grid(ONTOLOGY~., scales = 'free') # 以 ONTOLOGY 类型分开绘图
  print(p2)
  
  ## Tree
  # 注意：treeplot在某些版本中可能存在兼容性问题，使用tryCatch捕获错误
  tryCatch({
    GO_ed <- pairwise_termsim(GO)
    p3 <- treeplot(GO_ed)
    print(p3)
  }, error = function(e) {
    cat("警告：treeplot绘图失败，已跳过。错误信息：", e$message, "\n")
    cat("提示：这可能是包版本兼容性问题，dotplot和barplot已成功生成。\n")
  })
  dev.off()
  
}





# ----------------------------------- GO 结果筛选 ----------------------------------

# 1. 定义你想要保留的 GO ID 列表
# (你可以从你的 csv 文件中复制粘贴你感兴趣的 ID 到这里)
selected_IDs <- c(
  "GO:0042130", 
  "GO:0050868", 
  "GO:0032945"
  # ... 继续添加你想要的 ID
)

# 或者，如果你是根据关键词筛选的（比如只看 "T cell" 相关的）
# selected_rows <- grep("T cell", GO@result$Description)
# selected_IDs <- GO@result$ID[selected_rows]

# 2. 创建一个新的对象用于绘图（防止覆盖原始结果）
GO_filtered <- GO

# 3. 【核心步骤】将对象内部的数据替换为你筛选后的数据
# 这里的 @result 是 S4 对象存储数据的槽位
GO_filtered@result <- GO@result[GO@result$ID %in% selected_IDs, ]

# 4. 重新设置 pvalue 和 qvalue 的阈值检查
# 因为你已经手动筛选了，为了防止绘图函数内部再次过滤，我们可以把对象的阈值设宽
GO_filtered@pvalueCutoff <- 1
GO_filtered@qvalueCutoff <- 1

# 5. 检查筛选后还剩多少个
print(paste0("手动筛选后剩余条目数: ", nrow(GO_filtered@result)))

# 6. 使用筛选后的对象重新绘图
if (nrow(GO_filtered@result) > 0) {
  
  # 设置输出文件名
  pdf(file = paste0(name, "_GO_Manual_Filtered.pdf"), width = 10, height = 8)
  
  # --- 气泡图 ---
  # showCategory 可以设置得很大，以确保显示你选中的所有条目
  p1 <- dotplot(GO_filtered, 
                color = colorSel, 
                showCategory = nrow(GO_filtered@result), 
                title = "Manual Filtered GO Enrichment",
                label_format = 60,
                split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY~., scales = 'free')
  print(p1)
  
  # --- 条形图 ---
  p2 <- barplot(GO_filtered, 
                color = colorSel, 
                showCategory = nrow(GO_filtered@result), 
                title = "Manual Filtered GO Enrichment",
                label_format = 60,
                split = "ONTOLOGY") + 
        facet_grid(ONTOLOGY~., scales = 'free')
  print(p2)
  
  dev.off()
  print("筛选后的图片已生成。")
} else {
  print("错误：筛选后的 ID 在原始结果中找不到，请检查 ID 是否拼写正确。")
}





# ----------------------------------- clusterProfiler::simplify ----------------------------------


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(org.Rn.eg.db) # 大鼠，褐家鼠（Rat)

gene_list_all <- readRDS('./mfuzz_monocyte_6_gene.rds')

GO_draw_rat_gene_list <- function(gene_list_all){
  #参考：https://blog.csdn.net/qq_42090739/article/details/127306616
  #library(viridis)
  library(ggplot2)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(DOSE)
  library(cowplot)
  library(stringr)
  library(forcats)
  #org.Rn.eg.db <-loadDb("/hwfssz5/ST_SUPERCELLS/P21Z10200N0171/PROJECT/Skin_standby/01.ATAC_rn7/00.data/annotation/orgDB_Rattus.sqlite")
  
  # 提取基因和分组信息
  gene_name = unlist(gene_list_all)
  group_name = rep(names(gene_list_all), times = sapply(gene_list_all, length))
  df <- data.frame(gene = gene_name, group = group_name)
  
  data_GO <- compareCluster(
    gene~group, 
    data=df, 
    fun="enrichGO", 
    OrgDb= org.Rn.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    keyType = "SYMBOL"
  )
  
  # 富集分析结果简化:  表示仅保留相似性小于 0.7 的 GO 术语,根据 p.adjust p 值进行筛选
  data_GO_sim <- clusterProfiler::simplify(data_GO, 
                                           cutoff=0.7, 
                                           by="p.adjust", 
                                           select_fun=min)
  
  p = dotplot(data_GO_sim, showCategory=10, font.size = 8, label_format = 50)
  
  return(p)
  # ggsave("mono_all_GO.png", plot = p, width = 8, height = 15, dpi = 300)
  ggsave("mono_all_GO.pdf", plot = p, width = 8, height = 12)

}


gene_list_all <- readRDS('./mfuzz_monocyte_6_gene.rds')
gene_name = unlist(gene_list_all)
group_name = rep(names(gene_list_all), times = sapply(gene_list_all, length))
df <- data.frame(gene = gene_name, group = group_name)

data_GO <- compareCluster(
                          gene~group, 
                          data=df, 
                          fun="enrichGO", 
                          OrgDb= org.Rn.eg.db,
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          keyType = "SYMBOL")

# 富集分析结果简化:  表示仅保留相似性小于 0.7 的 GO 术语,根据 p.adjust p 值进行筛选
data_GO_sim <- clusterProfiler::simplify(data_GO, 
                                         cutoff=0.7, 
                                         by="p.adjust", 
                                         select_fun=min) # 多个相似术语中选择 p 值最小的术语

p = dotplot(data_GO_sim, showCategory=10, font.size = 8, label_format = 50)
# return(p)
ggsave("all_GO_simplify.pdf", plot = p, width = 8, height = 12)


# 提取和保存唯一的 GO 条目
unique_terms <- unique(data_GO_sim@compareClusterResult$Description)
write.csv(unique_terms, file = 'unique_term.csv', row.names = TRUE)

# 读取保存的条目并筛选 GO 结果
term <- read.csv('unique_term.csv', row.names = 1) 
filtered_GO <- data_GO_sim@compareClusterResult[data_GO_sim@compareClusterResult$Description %in% term$GO, ]

# 筛选出特定群集的 GO 条目
filtered_GO <- filtered_GO[filtered_GO$Cluster %in% c("cluster_1", "cluster_2", "cluster_3"), ]
data_GO_filtered <- data_GO_sim
data_GO_filtered@compareClusterResult <- as.data.frame(filtered_GO)
write.csv(data_GO_filtered@compareClusterResult, file = 'filtered_meta.csv', row.names = TRUE)

# 读取并重新排序处理后的数据
data_GO_filtered <- readRDS('data_GO_filtered.rds')
data_GO_filtered@compareClusterResult$Cluster <- factor(data_GO_filtered@compareClusterResult$Cluster,
                                                        levels = c("cluster_2", "cluster_1", "cluster_3"))

# 最终的绘图和保存
p <- dotplot(data_GO_filtered, showCategory = 10, font.size = 8, label_format = 50) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("selected_GO_filtered_END.pdf", plot = p, width = 5, height = 5.5)






# ----------------------------------- KEGG ----------------------------------

options(timeout = 300)
options(download.file.method = "libcurl")
options(url.method = "libcurl")


# 筛选上调表达的基因
top_up <- marker %>% filter(cluster == name) %>% 
          filter(avg_log2FC > 0.3 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) 
nrow(top_up)

# 筛选下调表达的基因
top_down <- marker %>% filter(cluster == name) %>% 
          filter(avg_log2FC < 0.3 & p_val_adj < 0.05) %>% arrange(avg_log2FC, p_val_adj) #%>% head(100)

nrow(top_down)

# 转换基因符号为ENTREZID
gene_ID_up <- bitr(top_up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
gene_ID_down <- bitr(top_down$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# KEGG富集分析
KEGG_up <- enrichKEGG(gene = gene_ID_up$ENTREZID, organism = "hsa", pvalueCutoff = pvalueFilter, qvalueCutoff = qvalueFilter)
print(paste0('The number of KEGG_up is : ',nrow(KEGG_up)))

KEGG_down <- enrichKEGG(gene = gene_ID_down$ENTREZID, organism = "hsa", pvalueCutoff = pvalueFilter, qvalueCutoff = qvalueFilter)
print(paste0('The number of KEGG_down is : ',nrow(KEGG_down)))



## Plots
if (nrow(KEGG_up) != 0) {
  
  pdf(file = paste0(name,"_KEGG_up.pdf"), width = 12, height = 10)
  
  p1 <- barplot(KEGG_up, x = "Count", color = colorSel, showCategory = 10, font.size = 15, title = "KEGG enrichment barplot", label_format = 50) # 超过40个字符串换行
  print(p1)
  p2 <- dotplot(KEGG_up, x = "GeneRatio", color = colorSel, showCategory = 10, title = "Top 10 of Pathway Enrichment",  label_format = 70) # 超过40个字符串换行
  print(p2)
  # 注意：treeplot在某些版本中可能存在兼容性问题，使用tryCatch捕获错误
  tryCatch({
    KEGG_up_ed <- pairwise_termsim(KEGG_up)
    p3 <- treeplot(KEGG_up_ed)
    print(p3)
  }, error = function(e) {
    cat("警告：KEGG_up treeplot绘图失败，已跳过。错误信息：", e$message, "\n")
    cat("提示：这可能是包版本兼容性问题，dotplot和barplot已成功生成。\n")
  })
  
  dev.off()
  
}


## Plots
if (nrow(KEGG_down) != 0 ) {
  
  pdf(file = paste0(name,"_KEGG_down.pdf"), width = 12, height = 10)
  
  p1 <- barplot(KEGG_down, x = "Count", color = colorSel, showCategory = 10, font.size = 15, title = "KEGG enrichment barplot", label_format = 50) # 超过40个字符串换行
  print(p1)
  p2 <- dotplot(KEGG_down, x = "GeneRatio", color = colorSel, showCategory = 10, title = "Top 10 of Pathway Enrichment",  label_format = 70) # 超过40个字符串换行
  print(p2)
  
    # 注意：treeplot在某些版本中可能存在兼容性问题，使用tryCatch捕获错误
    tryCatch({
      KEGG_up_ed <- pairwise_termsim(KEGG_down)
      p3 <- treeplot(KEGG_down_ed)
      print(p3)
    }, error = function(e) {
      cat("警告：KEGG_down treeplot 绘图失败，已跳过。错误信息：", e$message, "\n")
      cat("提示：这可能是包版本兼容性问题，dotplot和barplot已成功生成。\n")
    })
  
  dev.off()
  
}


# df = as.data.frame(KEGG_up)
# df$geneID = as.character(sapply(df$geneID, function(x) paste(gene_ID$SYMBOL[match(strsplit(x, "/")[[1]], as.character(gene_ID$ENTREZID))], collapse = "/"))) # 将geneID转换为基因符号
# df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
# write.csv(df_ed, file = paste0(name,"_KEGG_down.csv"), row.names = FALSE)
# 
# 
# df = as.data.frame(KEGG_down)
# df$geneID = as.character(sapply(df$geneID, function(x) paste(gene_ID$SYMBOL[match(strsplit(x, "/")[[1]], as.character(gene_ID$ENTREZID))], collapse = "/"))) # 将geneID转换为基因符号
# df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
# write.csv(df_ed, file = paste0(name,"_KEGG_down.csv"), row.names = FALSE)


###################### ------------------- 循环 ------------------- ######################
# # 开始循环处理每个细胞类型
# 开始循环处理每个细胞类型
# for (name in cell_types) {
  
#   # 打印当前处理状态
#   cat(paste0("Processing: ", name, "\n"))
  
#   # ------------------------ 基因筛选部分 ------------------------
#   tryCatch({
#     # 筛选上调表达的基因
#     top_up <- marker %>% 
#       filter(cluster == name) %>% 
#       filter(avg_log2FC > 2 & p_val < 0.05) %>% 
#       arrange(desc(scores), p_val_adj) %>% 
#       filter(!grepl("^ENSSS", gene) & 
#                !grepl("^ENSMFA", gene) & 
#                !grepl("^(RPL|RPS)", gene))
    
#     # 转换基因符号为 ENTREZID
#     gene_ID <- bitr(top_up$gene, 
#                     fromType = "SYMBOL", 
#                     toType = "ENTREZID", 
#                     OrgDb = org.Hs.eg.db)
    
#     # 如果没有有效基因则跳过
#     if(nrow(gene_ID) == 0) {
#       cat(paste0("No valid genes for ", name, "\n"))
#       next
#     }
    
#     # ------------------------ GO 富集分析 ------------------------
#     GO <- enrichGO(gene = gene_ID$ENTREZID,
#                    OrgDb = org.Hs.eg.db,
#                    keyType = "ENTREZID",
#                    ont = "ALL",
#                    pvalueCutoff = 0.5,
#                    qvalueCutoff = 0.5,
#                    readable = TRUE)
    
#     # 输出GO结果
#     # write.csv(as.data.frame(GO), file = paste0(name, "_GO_results.csv"))
    
#     # 绘制图形（如果有结果）
#     if (nrow(GO) != 0) {
#       # 生成GO图
#       pdf(file = paste0(name, "_GO.pdf"), width = 12, height = 15)
      
#       # 调整显示条目数
#       current_showNum <- ifelse(nrow(GO) < 10, nrow(GO), showNum)
      
#       p1 <- dotplot(GO, x = "GeneRatio", color = colorSel, size = "Count", 
#                     showCategory = current_showNum, label_format=150, split="ONTOLOGY") +
#         facet_grid(ONTOLOGY~., scales = 'free')
#       print(p1)
      
#       p2 <- barplot(GO, x = "Count", color = colorSel, 
#                     showCategory = current_showNum, label_format=150, split="ONTOLOGY") + 
#         facet_grid(ONTOLOGY~., scales = 'free')
#       print(p2)
      
#       GO_ed <- pairwise_termsim(GO)
#       p3 <- treeplot(GO_ed)
#       print(p3)
      
#       dev.off()
#     }
    
#     # ------------------------ KEGG 富集分析 ------------------------
#     KEGG_up <- enrichKEGG(gene = gene_ID$ENTREZID, 
#                           organism = "hsa", 
#                           pvalueCutoff = pvalueFilter, 
#                           qvalueCutoff = qvalueFilter)
    
#     # 输出KEGG结果
#     # write.csv(as.data.frame(KEGG_up), file = paste0(name, "_KEGG_results.csv"))
    
#     # 绘制图形（如果有结果）
#     if (nrow(KEGG_up) != 0) {
#       # 生成KEGG图
#       pdf(file = paste0(name, "_KEGG_up.pdf"), width = 12, height = 10)
      
#       # 调整显示条目数
#       current_showNum <- ifelse(nrow(KEGG_up) < 10, nrow(KEGG_up), showNum)
      
#       p1 <- barplot(KEGG_up, x = "Count", color = colorSel, 
#                     showCategory = current_showNum, font.size = 15, 
#                     title = "KEGG enrichment barplot", label_format = 50)
#       print(p1)
      
#       p2 <- dotplot(KEGG_up, x = "GeneRatio", color = colorSel, 
#                     showCategory = current_showNum, 
#                     title = paste0("Top ", current_showNum ," of Pathway Enrichment"),  
#                     label_format = 70)
#       print(p2)
      
#       KEGG_up_ed <- pairwise_termsim(KEGG_up)
#       p3 <- treeplot(KEGG_up_ed)
#       print(p3)
      
#       dev.off()
#     }
    
#   }, error = function(e) {
#     cat(paste0("Error processing ", name, ": ", e$message, "\n"))
#   })
# }


# ----------------------------------- GSEA ----------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)
library(Seurat)
library(dplyr)

marker_gsea <- read.csv("./LUSC-Basal-Cell_markers.csv")
colnames(marker_gsea)[1] <- "SYMBOL" # 将第一列命名为'gene'

# 转换基因符号为ENTREZID并排序
genelist <- marker_gsea %>% 
  select(SYMBOL, avg_log2FC) %>%
  merge(bitr(.$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db), by = "SYMBOL") %>%
  arrange(desc(avg_log2FC)) %>%
  {setNames(.$avg_log2FC, .$ENTREZID)}
head(genelist)


# GO富集分析
gsego <- gseGO(geneList     = genelist,#根据LogFC排序后的基因列表
               OrgDb        = org.Hs.eg.db,
               ont          = "ALL",# GO分析的模块
               minGSSize    = 10, #最小基因集的基因数
               maxGSSize    = 500,#最大基因集的基因数
               pvalueCutoff = 0.05,#p值的阈值
               verbose      = FALSE)#是否输出提示信息

df_go <- as.data.frame(gsego)
df_go_filtered <- df_go[df_go$pvalue < pvalueFilter & df_go$qvalue < qvalueFilter, ]
write.csv(df_go_filtered, file = paste0(name, "_gseGO.csv"), row.names = FALSE)

p1 <- gseaplot2(gsego, geneSetID = 1:3)
p1
ggsave("gseGO.png", p1, bg = "white", dpi = 300, width = 8, height = 7)


# KEGG富集分析
gseKEGG <- gseKEGG(geneList     = genelist,#根据LogFC排序后的基因列表
                   organism = "hsa",
                   pvalueCutoff = 0.05,# p值的阈值
                   verbose      = FALSE)#是否输出提示信息

df_kegg <- as.data.frame(gseKEGG)
df_kegg_filtered <- df_kegg[df_kegg$pvalue < pvalueFilter & df_kegg$qvalue < qvalueFilter, ]
write.csv(df_kegg_filtered, file = paste0(name, "_gseKEGG.csv"), row.names = FALSE)

plot_gseKEGG <- gseaplot2(gseKEGG, geneSetID = c("hsa04657","hsa05323"))
plot_gseKEGG
ggsave("gseKEGG.png", plot_gseKEGG, bg = "white", dpi = 300, width = 8, height = 7)


























