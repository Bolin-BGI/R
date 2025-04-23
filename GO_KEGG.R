library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(GOplot)



# 读取CSV文件，去除 'p-' 'm-' 前缀
marker <- read.csv("./YWHAG_markers.csv")
colnames(marker)[1] <- "gene" # 将第一列命名为'gene'


# 筛选上调表达的基因并选择前 100个基因
top_up <- marker %>% filter(avg_log2FC > 0.3 & p_val < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) %>% 
          filter(!grepl("^ENSSS", gene) & !grepl("^ENSMFA", gene) & !grepl("^(RPL|RPS)", gene)) # %>% head(100)

# 转换基因符号为ENTREZID
gene_ID <- bitr(top_up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# ----------------------------------- GO ----------------------------------

# 设置p值和q值过滤阈值
pvalueFilter = 0.05  
qvalueFilter = 1

# 根据q值过滤阈值选择颜色
colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

# main
GO <- enrichGO(gene = gene_ID$ENTREZID,OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",  # 记得和输入基因的格式保持一致
              ont = "ALL", # BP(Biological Process)  MF(Molecular Function)  CC(Cellular Component) 
              pvalueCutoff = 0.5, 
              qvalueCutoff = 0.5,
              readable = TRUE)


df = as.data.frame(GO)
df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
write.csv(df_ed, file = "YWHAG_positive_GO_up.csv", row.names = FALSE)

# 设置显示的条目数
showNum = 10
if (nrow(GO) < 10) {
  showNum = nrow(GO)
}

# 如果有富集结果，生成气泡图
if (nrow(GO) != 0) {
  # 生成气泡图
  pdf(file = paste0(name, "_GO.pdf"), width = 10, height = 15)
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

# 设置p值和q值过滤阈值
pvalueFilter = 0.2 # 0.05  
qvalueFilter = 0.2 # 1

# 根据q值过滤阈值选择颜色
colorSel = "qvalue"
if (qvalueFilter > 0.05) {
  colorSel = "pvalue"
}

# 筛选上调表达的基因并选择前 100个基因
top_up <- marker %>% filter(avg_log2FC > 0.3 & p_val_adj < 0.05) %>% arrange(desc(avg_log2FC), p_val_adj) 
# 筛选下调表达的基因并选择前 100 个基因
top_down <- marker %>% filter(avg_log2FC < 0.3 & p_val_adj < 0.05) %>% arrange(avg_log2FC, p_val_adj) #%>% head(100)

# 转换基因符号为ENTREZID
gene_ID_up <- bitr(top_up$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
gene_ID_down <- bitr(top_down$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG富集分析
KEGG_up <- enrichKEGG(gene = gene_ID_up$ENTREZID, organism = "hsa", pvalueCutoff = pvalueFilter, qvalueCutoff = qvalueFilter)
KEGG_down <- enrichKEGG(gene = gene_ID_down$ENTREZID, organism = "hsa", pvalueCutoff = pvalueFilter, qvalueCutoff = qvalueFilter)



## Plots
if (nrow(KEGG_up) != 0) {
  
  pdf(file = paste0(name,"_KEGG_up.pdf"), width = 10, height = 10)
  
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
if (nrow(KEGG_down) != 0) {
  
  pdf(file = paste0(name,"_KEGG_down.pdf"), width = 10, height = 10)
  
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



df = as.data.frame(KEGG_up)
df$geneID = as.character(sapply(df$geneID, function(x) paste(gene_ID$SYMBOL[match(strsplit(x, "/")[[1]], as.character(gene_ID$ENTREZID))], collapse = "/"))) # 将geneID转换为基因符号
df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
write.csv(df_ed, file = "YWHAG_KEGG_up.csv", row.names = FALSE)


df = as.data.frame(KEGG_down)
df$geneID = as.character(sapply(df$geneID, function(x) paste(gene_ID$SYMBOL[match(strsplit(x, "/")[[1]], as.character(gene_ID$ENTREZID))], collapse = "/"))) # 将geneID转换为基因符号
df_ed = df[(df$pvalue < pvalueFilter & df$qvalue < qvalueFilter),] 
write.csv(df_ed, file = "YWHAG_KEGG_down.csv", row.names = FALSE)






