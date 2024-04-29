B/Plasma Figure
================

## Set up

Load R libraries

``` r
library(ggplot2)
library(ggpubr)
library(glue)
library(gtools)
library(tidyverse)

library(reticulate)
use_python("/projects/home/tlchan/.conda/envs/myenv/bin/python")
```

Load python libraries

``` python
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pegasus as pg
import scanpy as sc
```

## Figure 1A

``` python
mpl.rcParams['pdf.fonttype'] = 42

cluster_palette = {
    "1": "#FF0029",
    "2": "#377EB8",
    "3": "#66A61E",
    "4": "#984EA3",
    "5": "#00D2D5",
    "6": "#FF7F00",
    "7": "#AF8D00",
    "8": "#7F80CD",
    "9": "#B3E900",
    "10": "#C42E60"
}

# Load single-cell object
lineage_data = pg.read_input(
    '/projects/home/tlchan/projects/ascites/second_data_freeze/clusterings/ascites_bplasma_R5_300mg_20pm_harm_channel_multi_res/1.5/data/pseudobulk/ascites_bplasma_R5_300mg_20pm_harm_channel_1_5_complete_with_pb.zarr.zip')

# Relabel obs for function
lineage_data.obs['Cluster'] = lineage_data.obs['leiden_labels'].cat.remove_unused_categories().astype(str)

cluster_fig, cluster_ax = plt.subplots(1)
cluster_umap = sc.pl.umap(adata=lineage_data.to_anndata(),
                          color="Cluster",
                          use_raw=True,
                          palette=cluster_palette,
                          legend_loc="on data",
                          legend_fontoutline=5,
                          title="",
                          show=False,
                          ax=cluster_ax)

cluster_fig = plt.gcf()
cluster_fig.set_size_inches(6, 6)
cluster_fig.tight_layout()
cluster_ax.set_rasterization_zorder(2)
plt.show()
plt.close()
```

    ## 2024-04-25 17:33:43,361 - pegasusio.readwrite - INFO - zarr file '/projects/home/tlchan/projects/ascites/second_data_freeze/clusterings/ascites_bplasma_R5_300mg_20pm_harm_channel_multi_res/1.5/data/pseudobulk/ascites_bplasma_R5_300mg_20pm_harm_channel_1_5_complete_with_pb.zarr.zip' is loaded.
    ## 2024-04-25 17:33:43,362 - pegasusio.readwrite - INFO - Function 'read_input' finished in 0.38s.
    ## /projects/home/tlchan/.conda/envs/myenv/lib/python3.9/site-packages/scanpy/plotting/_tools/scatterplots.py:392: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
    ##   cax = scatter(

<img src="bplasma_figure_files/figure-gfm/fig_1A-1.png" width="576" />

## Figure 1B

``` r
b_genes <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_markers.csv') %>%
    filter(lineage == 'bplasma') %>%
    pull(genes) %>%
    strsplit(",") %>%
    unlist()
b_gex <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_data/bplasma_gene_exp.csv') %>%
    mutate(Gene = factor(Gene, levels = b_genes), Cluster = factor(Cluster))

b_gp <- ggplot(b_gex, aes(x = Gene, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
    geom_point(color = "black", shape = 21) +
    ylab('Cluster') +
    labs(size = "% Expressed") +
    scale_fill_gradient(low = "#fff5f0", high = "#67000c") +
    lims(size = c(0, 100)) +
    theme_light(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(legend.key.size = unit(.25, "cm"))

b_proteins <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_markers.csv') %>%
    filter(lineage == 'bplasma') %>%
    pull(proteins) %>%
    strsplit(",") %>%
    unlist()

b_cite <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/dotplot_data/bplasma_cite_exp.csv') %>%
    mutate(Protein = factor(Protein, levels = b_proteins), Cluster = factor(Cluster))

b_pp <- ggplot(b_cite, aes(x = Protein, y = fct_rev(Cluster), fill = Count, size = Percent_Expressed)) +
    geom_point(color = "black", shape = 21) +
    labs(size = "% Expressed") +
    scale_fill_gradient(low = "#fff5f0", high = "#08306b") +
    lims(size = c(0, 100)) +
    theme_light(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    theme(legend.key.size = unit(.25, "cm"))

ggarrange(b_gp, b_pp, ncol = 2, nrow = 1, widths = c(1.0, 0.4), align = "h")
```

![](/tmp/bplasma_figure-2.rmd/bplasma_figure_files/figure-gfm/fig_1B-3.png)<!-- -->
