library(ggpubr)
library(tidyverse)

# All dotplot markers
dotplot_markers <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_markers.csv')

plot_dotplot <- function(lin_gex, lin_cite, lin, widths) {
    lin_genes <- dotplot_markers %>%
        filter(lineage == lin) %>%
        pull(genes) %>%
        strsplit(",") %>%
        unlist()
    lin_gex <- lin_gex %>%
        mutate(Gene = factor(Gene, levels = lin_genes), Cluster = factor(Cluster))

    lin_gp <- ggplot(lin_gex, aes(x = Gene, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
        geom_point(color = "black", shape = 21) +
        ylab('Cluster') +
        labs(size = "% Expressed") +
        scale_fill_gradient(low = "#fff5f0", high = "#67000c") +
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

    lin_pp <- ggplot(lin_cite, aes(x = Protein, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
        geom_point(color = "black", shape = 21) +
        labs(size = "% Expressed") +
        scale_fill_gradient(low = "#fff5f0", high = "#08306b") +
        lims(size = c(0, 100)) +
        theme_light(base_size = 25) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
        theme(legend.key.size = unit(.4, "cm"), legend.title = element_text(size = 20), legend.text = element_text(size = 12))

    legends <- ggarrange(get_legend(lin_gp), get_legend(lin_pp), nrow = 2, align = "v")

    ggarrange(lin_gp + theme(legend.position = "none"), lin_pp + theme(legend.position = "none"), legends, ncol = 3, nrow = 1, widths = widths, align = "h")
}