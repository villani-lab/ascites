devtools::install("/projects/home/tlchan/util/UpSetR", force = TRUE)
library(tidyverse)
library(UpSetR)

plot_upset <- function(res) {
    # Select upregulated data
    sig_up <- res[(res$padj < 0.1) & (res$log2FoldChange > 0),]

    # Lets make a list for the upset plot
    upset_list_up <- sig_up %>%
        split(.$cluster) %>%
        lapply(., function(x) x$gene_symbol)

    upset(fromList(upset_list_up),
          sets = names(upset_list_up),
          scale.intersections = "log10",
          order.by = "freq",
          point.size = 3,
          line.size = 1.5,
          keep.order = TRUE,
          sets.x.label = "# genes",
          mainbar.y.label = "# overlapping genes",
          text.scale = c(1.6, 1.8, 1.6, 1.8, 1.8, 0))
}