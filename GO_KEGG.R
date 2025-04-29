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
if (nrow(GO) < 10) {
  showNum = nrow(GO)
}

# main
GO <- enrichGO(gene = gene_ID$ENTREZID,OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",  # 记得和输入基因的格式保持一致
               ont = "ALL", # BP(Biological Process)  MF(Molecular Function)  CC(Cellular Component) 
               pvalueCutoff = 0.5, 
               qvalueCutoff = 0.5,
               readable = TRUE)
print(paste0('The number of GO_items is : ',nrow(GO)))

# df = as.data.frame(GO)
# df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
# write.csv(df_ed, file = "YWHAG_positive_GO_up.csv", row.names = FALSE)

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
  GO_ed <- pairwise_termsim(GO)
  p3 <- treeplot(GO_ed)
  print(p3)
  dev.off()
  
}



# ----------------------------------- KEGG ----------------------------------

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
  ## Tree
  KEGG_up_ed <- pairwise_termsim(KEGG_up)
  p3 <- treeplot(KEGG_up_ed)
  print(p3)
  
  dev.off()
  
}


## Plots
if (nrow(KEGG_down) != 0 ) {
  
  pdf(file = paste0(name,"_KEGG_down.pdf"), width = 12, height = 10)
  
  p1 <- barplot(KEGG_down, x = "Count", color = colorSel, showCategory = 10, font.size = 15, title = "KEGG enrichment barplot", label_format = 50) # 超过40个字符串换行
  print(p1)
  p2 <- dotplot(KEGG_down, x = "GeneRatio", color = colorSel, showCategory = 10, title = "Top 10 of Pathway Enrichment",  label_format = 70) # 超过40个字符串换行
  print(p2)
  ## Tree
  KEGG_down_ed <- pairwise_termsim(KEGG_down)
  p3 <- treeplot(KEGG_down_ed)
  print(p3)
  
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
#     gene_ID <- bitr(top_up$gene, fromType = "SYMBOL", toType = "ENTREZID",  OrgDb = org.Hs.eg.db)
    
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
#     if(nrow(GO) > 0) {
#       current_showNum <- min(showNum, nrow(GO))
      
#       pdf(file = paste0(name, "_GO_dotplot.pdf"), width = 12, height = 15)
#       dotplot(GO, x = "GeneRatio", color = colorSel, size = "Count", 
#               showCategory = current_showNum, label_format=150, split="ONTOLOGY") + 
#         facet_grid(ONTOLOGY~., scales = 'free')
#       dev.off()
      
#       pdf(file = paste0(name, "_GO_barplot.pdf"), width = 12, height = 8)
#       barplot(GO, x = "Count", color = colorSel, 
#               showCategory = current_showNum, label_format=150, split="ONTOLOGY") + 
#         facet_grid(ONTOLOGY~., scales = 'free')
#       dev.off()
      
#       # 树状图
#       GO_ed <- pairwise_termsim(GO)
#       pdf(file = paste0(name, "_GO_treeplot.pdf"))
#       treeplot(GO_ed)
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
#     if(nrow(KEGG_up) > 0) {
#       current_showNum <- min(showNum, nrow(KEGG_up))
      
#       pdf(file = paste0(name, "_KEGG_barplot.pdf"), width = 12, height = 10)
#       barplot(KEGG_up, x = "Count", color = colorSel, 
#               showCategory = current_showNum, font.size = 15, 
#               title = "KEGG enrichment barplot", label_format = 50)
#       dev.off()
      
#       pdf(file = paste0(name, "_KEGG_dotplot.pdf"), width = 12, height = 10)
#       dotplot(KEGG_up, x = "GeneRatio", color = colorSel, 
#               showCategory = current_showNum, 
#               title = paste0("Top ", current_showNum ," of Pathway Enrichment"),  
#               label_format = 70)
#       dev.off()
      
#       # 树状图
#       KEGG_up_ed <- pairwise_termsim(KEGG_up)
#       pdf(file = paste0(name, "_KEGG_treeplot.pdf"))
#       treeplot(KEGG_up_ed)
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


























