library(ggpubr)
library(tidyverse)

# Use Reds
reds_colormap <- read.table('/projects/home/tlchan/util/cmaps/reds_cmap.txt', header = FALSE)
reds_hex <- apply(reds_colormap, 1, function(row) {
    rgb(row[1], row[2], row[3], maxColorValue = 1)
})

# Use Blues
blues_colormap <- read.table('/projects/home/tlchan/util/cmaps/blues_cmap.txt', header = FALSE)
blues_hex <- apply(blues_colormap, 1, function(row) {
    rgb(row[1], row[2], row[3], maxColorValue = 1)
})

# All dotplot markers
dotplot_markers <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_markers.csv')

plot_dotplot <- function(lin_gex, lin_cite, lin, widths, cluster_order) {
    lin_genes <- dotplot_markers %>%
        filter(lineage == lin) %>%
        pull(genes) %>%
        strsplit(",") %>%
        unlist()
    lin_gex <- lin_gex %>%
        mutate(Gene = factor(Gene, levels = lin_genes), Cluster = factor(Cluster))

    if (hasArg(cluster_order)) {
        lin_gex <- lin_gex %>% mutate(Cluster = factor(Cluster, levels = cluster_order))
    }

    lin_gp <- ggplot(lin_gex, aes(x = Gene, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
        geom_point(color = "black", shape = 21) +
        ylab('Cluster') +
        labs(size = "% Expressed", fill = "Scaled expression") +
        scale_fill_gradientn(colors = reds_hex) +
        lims(size = c(0, 100)) +
        theme_light(base_size = 25) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(legend.key.size = unit(.4, "cm"), legend.title = element_text(size = 20), legend.text = element_text(size = 12)) +
        guides(size = "none")

    lin_proteins <- dotplot_markers %>%
        filter(lineage == lin) %>%
        pull(proteins) %>%
        strsplit(",") %>%
        unlist()

    lin_cite <- lin_cite %>%
        mutate(Protein = factor(Protein, levels = lin_proteins), Cluster = factor(Cluster))

    if (hasArg(cluster_order)) {
        lin_cite <- lin_cite %>% mutate(Cluster = factor(Cluster, levels = cluster_order))
    }

    lin_pp <- ggplot(lin_cite, aes(x = Protein, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
        geom_point(color = "black", shape = 21) +
        labs(size = "% Expressed", fill = "Scaled expression") +
        scale_fill_gradientn(colors = blues_hex) +
        lims(size = c(0, 100)) +
        theme_light(base_size = 25) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        theme(legend.key.size = unit(.4, "cm"), legend.title = element_text(size = 20), legend.text = element_text(size = 12))

    legends <- ggarrange(get_legend(lin_gp), get_legend(lin_pp), nrow = 2, align = "v")

    ggarrange(lin_gp + theme(legend.position = "none"), lin_pp + theme(legend.position = "none"), legends, ncol = 3, nrow = 1, widths = widths, align = "h")
}