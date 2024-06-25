Supplemental CD4 Figure
================

## Set up

Load R libraries

``` r
library(tidyverse)

library(reticulate)
use_python("/projects/home/tlchan/.conda/envs/myenv/bin/python")

setwd('/projects/home/tlchan/github_code/ascites/functions')
source('plot_upset.R')
```

Load python libraries

``` python
import pegasus as pg

import sys
sys.path.append("/projects/home/tlchan/github_code/ascites/functions")
import python_functions
```

## Figure 1A

``` python
cd4_data = pg.read_input("/projects/home/tlchan/projects/ascites/data_cite_objects/cd4.zarr.zip")

python_functions.plot_feature(lin_data=cd4_data,
                              genes=['GZMB', 'GZMK', 'TBX21', 'cite_CD45RA', 'IFNG', 'CD8A', 'LAIR2', 'cite_CD45RO'],
                              ncol=4,
                              nrow=2)
```

    ## 2024-06-20 18:49:32,980 - pegasusio.readwrite - INFO - zarr file '/projects/home/tlchan/projects/ascites/data_cite_objects/cd4.zarr.zip' is loaded.
    ## 2024-06-20 18:49:32,980 - pegasusio.readwrite - INFO - Function 'read_input' finished in 0.84s.

<img src="supp_cd4_figure_files/figure-gfm/fig_1A-1.png" width="1920" />

## Figure 1B

``` r
cd4_cluster_palette <- list("1" = "#FF0029",
                            "2" = "#377EB8",
                            "3" = "#66A61E",
                            "4" = "#984EA3",
                            "5" = "#00D2D5",
                            "6" = "#FF7F00",
                            "7" = "#AF8D00",
                            "8" = "#7F80CD",
                            "9" = "#B3E900"
)

cd4_data <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/abundance_data/cd4_patient_counts.csv') %>%
    mutate(Cluster = factor(Cluster))

ggplot(cd4_data, aes(x = Patient, y = Count, fill = Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = cd4_cluster_palette)
```

![](/tmp/supp_cd4_figure-3.rmd/supp_cd4_figure_files/figure-gfm/fig_1B-3.png)<!-- -->

## Figure 1C

``` r
cd4_res <- read.csv('/projects/home/tlchan/projects/ascites/ascites_deg_results/cd4/cd4_de_by_tissue_type_all_results.csv')

plot_upset(res = cd4_res,
           min_degree = 7)
```

![](/tmp/supp_cd4_figure-3.rmd/supp_cd4_figure_files/figure-gfm/fig_1C-1.png)<!-- -->
