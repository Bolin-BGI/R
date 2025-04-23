colorlist = c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941","#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC","#B79762", "#004D43", "#8FB0FF", 
              "#997D87", "#5A0007","#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80","#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9","#B903AA", 
              "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8","#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400","#A1CAF1", "#C2B280", "#848482", 
              "#E68FAC", "#0067A5","#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300","#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26","#191970", 
              "#000080","#6495ED", "#1E90FF", "#00BFFF", "#00FFFF", "#FF1493","#FF00FF", "#A020F0", "#63B8FF", "#008B8B", "#54FF9F","#00FF00", "#76EE00", "#FFF68F")

# if (length(Cells(combined))>60000) {
#     draw_obj = subset(combined, cells=sample(Cells(combined), 50000))
# } else {
#     draw_obj = combined
# }

plots <- function(obj, group_name, prefix, ...){
    
    colors = colorlist[1:length(unique(obj@meta.data[[group_name]]))]
    names(colors) = unique(obj@meta.data[[group_name]])
    
    p <- DimPlot(obj, group.by=group_name, pt.size=0.5, label.size = 8, cols=colors, raster=F, ...)
    ggsave(paste0(prefix, "_", group_name, ".png"),width=15,height=14, dpi = 300)
    
}

plots_split <- function(obj, draw_obj, group_name, prefix, ...){
    
    if (length(Cells(obj))>60000) {
    draw_obj = subset(obj, cells=sample(Cells(obj), 50000))
    } else {
        draw_obj = obj
    }
    
    colors = colorlist[1:length(unique(obj@meta.data[[group_name]]))]
    names(colors) = unique(obj@meta.data[[group_name]])
    
    pdf(paste0(prefix, "_", group_name, "_split.pdf"))
    for (i in sort(unique(draw_obj@meta.data[[group_name]]))){
        cells = Cells(draw_obj)[draw_obj@meta.data[[group_name]]==i]
        p <- DimPlot(draw_obj, group.by=group_name, cells.highlight=cells)
        print(p+ggtitle(i))
    }
    dev.off()
}
