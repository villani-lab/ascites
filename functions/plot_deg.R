library(ggrepel)
library(tidyverse)

plot_deg <- function(res, meta, lin, clust, contrast = "", ref_var = "", test_var = "", cell_cutoff = 5, logfc = 0, padj = 0.1) {
    meta_cluster <- meta[meta$cluster == clust,]

    # Filtering for sample/clusters that have more than cell_cutoff cells
    meta_cluster <- meta_cluster[meta_cluster$n_cells >= cell_cutoff,]
    tot_cells <- sum(meta_cluster$n_cells)

    # Filter results
    res <- res %>% filter(cluster == clust)

    # Make plot
    top <- res[res$log2FoldChange > logfc & res$padj < padj,]
    top20 <- head(top[order(top$padj),], 20L)
    bottom <- res[res$log2FoldChange < -logfc & res$padj < padj,]
    bottom20 <- head(bottom[order(bottom$padj),], 20L)
    if (contrast != "") {
        plot_title <- sprintf('%s %s: %s (%i) vs %s (%i)', lin, clust, test_var,
                              nrow(meta[meta[contrast] == test_var,]), ref_var,
                              nrow(meta[meta[contrast] == ref_var,]))
        plot_subtitle <- sprintf('%i total cells, abs(log2fc) > %1.1f, and padj < %1.2f',
                                 tot_cells, logfc, padj)
    } else {
        plot_title <- sprintf('%s %s', lin, clust)
        plot_subtitle <- sprintf('%i total cells, abs(log2fc) > %1.1f, and padj < %1.2f',
                                 tot_cells, logfc, padj)
    }

    p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(data = res[abs(res$log2FoldChange) < logfc | res$padj > padj,], color = "grey") +
        geom_point(data = top, color = "red") +
        geom_point(data = bottom, color = "blue") +
        geom_hline(yintercept = -log10(padj), linetype = 'dashed', size = 0.75) +
        geom_vline(xintercept = logfc, linetype = 'dashed', size = 0.75) +
        geom_vline(xintercept = -logfc, linetype = 'dashed', size = 0.75) +
        geom_text_repel(data = top20, aes(label = gene_symbol), max.overlaps = Inf, size = 5) +
        geom_text_repel(data = bottom20, aes(label = gene_symbol), max.overlaps = Inf, size = 5) +
        xlab("log2 Fold Change") +
        ggtitle(plot_title, subtitle = plot_subtitle) +
        theme_bw(base_size = 15)

    return(p)
}
