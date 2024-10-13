library(ggrepel)
library(tidyverse)

plot_lineage_deg <- function(res, meta, lin, contrast, ref_var, test_var, cell_cutoff = 10, fc_cutoff = 0, p_cutoff = 0.1, ref_color = "#2278B5", test_color = "#2278B5") {
    # Filtering for sample/clusters that have more than cell_cutoff cells
    meta <- meta[meta$n_cells >= cell_cutoff,]
    tot_cells <- sum(meta$n_cells)

    # Make plot
    top <- res[res$log2FoldChange > fc_cutoff & res$padj < p_cutoff,]
    top20 <- head(top[order(top$padj),], 20L)
    bottom <- res[res$log2FoldChange < -fc_cutoff & res$padj < p_cutoff,]
    bottom20 <- head(bottom[order(bottom$padj),], 20L)

    if (hasArg(contrast)) {
        plot_title <- sprintf('%s: %s (%i) vs %s (%i)', lin, test_var,
                              nrow(meta[meta[contrast] == test_var,]), ref_var,
                              nrow(meta[meta[contrast] == ref_var,]))
        plot_subtitle <- sprintf('%i total cells, %i cell cutoff, abs(log2fc) > %1.1f, and padj < %1.2f',
                                 tot_cells, cell_cutoff, fc_cutoff, p_cutoff)
    } else {
        plot_title <- sprintf('%s', lin)
        plot_subtitle <- sprintf('%i total cells, %i cell cutoff, abs(log2fc) > %1.1f, and padj < %1.2f',
                                 tot_cells, cell_cutoff, fc_cutoff, p_cutoff)
    }

    ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(data = res[abs(res$log2FoldChange) < fc_cutoff | res$padj > p_cutoff,], color = "grey") +
        geom_point(data = top, color = test_color) +
        geom_point(data = bottom, color = ref_color) +
        geom_hline(yintercept = -log10(p_cutoff), linetype = 'dashed', size = 0.75) +
        geom_vline(xintercept = fc_cutoff, linetype = 'dashed', size = 0.75) +
        geom_vline(xintercept = -fc_cutoff, linetype = 'dashed', size = 0.75) +
        geom_text_repel(data = top20, aes(label = gene_symbol), max.overlaps = Inf, size = 5) +
        geom_text_repel(data = bottom20, aes(label = gene_symbol), max.overlaps = Inf, size = 5) +
        xlab("log2 Fold Change") +
        ggtitle(plot_title, subtitle = plot_subtitle) +
        theme_bw(base_size = 15)
}