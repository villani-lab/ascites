library(tidyverse)
library(ComplexUpset)

plot_upset <- function(res, min_degree) {
    # Select upregulated data
    sig_up <- res %>% filter(padj < 0.05, log2FoldChange > 0)

    # Get list of significant genes per cluster
    upset_list_up <- sig_up %>%
        split(.$cluster) %>%
        lapply(., function(x) x$gene_symbol)


    # Prepare dataframe for upset plot
    clusters <- names(upset_list_up)
    gene_list <- unique(unlist(upset_list_up))
    genes <- unlist(lapply(upset_list_up, function(x) { x <- as.vector(match(gene_list, x)) }))
    genes[is.na(genes)] <- as.integer(0); genes[genes != 0] <- as.integer(1)
    genes <- matrix(genes, ncol = length(upset_list_up), byrow = F)
    genes <- data.frame(genes == 1)
    names(genes) <- clusters

    upset(genes,
          colnames(genes),
          name = 'Clusters',
          base_annotations = list(
              'Intersection size' = intersection_size(counts = FALSE) +
                  ylab("# overlapping genes") +
                  theme(text = element_text(size = 22))
          ),
          set_sizes = (
              upset_set_size() +
                  ylab("# genes") +
                  scale_y_reverse(n.breaks = 3) +
                  theme(text = element_text(size = 22))
          ),
          themes = upset_default_themes(text = element_text(size = 22)),
          width_ratio = 0.25,
          height_ratio = 1,
          min_degree = min_degree,
          sort_sets = FALSE)
}