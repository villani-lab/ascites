Supplemental Figure 3
================

## Set up

Load R libraries

``` r
library(reticulate)
use_python("/projects/home/tlchan/.conda/envs/myenv/bin/python")

setwd('/projects/home/tlchan/github_code/ascites/functions')
source('plot_abundance.R')
source('plot_dotplot.R')
```

Load python libraries

``` python
import matplotlib.pyplot as plt
import pegasus as pg

import sys
sys.path.append("/projects/home/tlchan/github_code/ascites/functions")
import python_functions
```

## Supplementary Figure 3A

``` python
b_cluster_palette = {
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
b_data = pg.read_input(
    '/projects/home/tlchan/projects/ascites/integrated_data/clusterings/ascites_bplasma_R5_300mg_20pm_harm_channel_multi_res/1.5/data/pseudobulk/ascites_bplasma_R5_300mg_20pm_harm_channel_1_5_complete_with_pb.zarr.zip')

# Relabel obs for function
b_data.obs['Cluster'] = b_data.obs['leiden_labels'].cat.remove_unused_categories().astype(str)

fig = python_functions.plot_umap(lin_data=b_data,
                           palette=b_cluster_palette)

plt.show()
# plt.savefig("/projects/home/tlchan/fig_panels/supp_3a.pdf")
plt.close(fig)
```

<img src="supp_figure_3_files/figure-gfm/supp_3A-1.png" width="576" />

## Supplementary Figure 3B

``` r
b_gex <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/data/dotplot_data/bplasma_gene_exp.csv')
b_cite <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/data/dotplot_data/bplasma_cite_exp.csv')

plot_dotplot(lin_gex = b_gex,
             lin_cite = b_cite,
             lin = "bplasma",
            widths=c(1, .15, .2))
```

![](supp_figure_3_files/figure-gfm/supp_3B-3.png)<!-- -->

``` r
# ggsave("/projects/home/tlchan/fig_panels/supp_3b.pdf", width = 12, height = 8)
```

## Supplementary Figure 3C

``` r
plot_cluster_abundance("bplasma")
```

![](supp_figure_3_files/figure-gfm/supp_3C-1.png)<!-- -->

``` r
# ggsave("/projects/home/tlchan/fig_panels/supp_3c.pdf", width = 8, height = 8)
```
