library(pheatmap)
library(RColorBrewer)

draw_heatmap = function(data1, data2, feature_file){
    mat1 = read.table(data1, header=TRUE, row.names=1)
    mat2 = read.table(data2, header=TRUE, row.names=1)

    locus = rownames(mat1)
    
    feature = read.table(feature_file, header=TRUE, row.names=1)[locus,]
    gc_sort_index = rownames(feature)[sort(feature$GC, index.return=TRUE)$ix]

    mat = cbind(mat1[gc_sort_index,], mat2[gc_sort_index,])/feature[gc_sort_index,"Length"]

    p = pheatmap(mat, color=rev(colorRampPalette(brewer.pal(n=9, "YlGnBu"))(100)),
                 show_colnames=TRUE, show_rownames=FALSE, cluster_cols=TRUE, cluster_rows=FALSE)
    return(p)
}

draw_heatmap("../ice_4_hapmap_exomes.covered_bases.txt", "../twist_2_hapmap_exomes.covered_bases.txt", "../gc_length_metadata.txt")
pheatmap(data.frame(feature[gc_sort_index, "GC"]), cluster_rows=F, cluster_cols=F, color=rev(heat.colors(100)), show_rownames=FALSE, show_colnames=FALSE)
