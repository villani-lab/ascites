Supplemental CD8 Figure
================

## Set up

Load R libraries

``` r
library(tidyverse)

library(reticulate)
use_python("/projects/home/tlchan/.conda/envs/myenv/bin/python")
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
cd8_data = pg.read_input("/projects/home/tlchan/projects/ascites/data_cite_objects/cd8.zarr.zip")

python_functions.plot_feature(lin_data=cd8_data,
                              genes=['cite_CD45RA', 'cite_CLEC12A', 'CLEC12A', 'cite_CD45RO', 'cite_CD11c', 'ITGAX'],
                              ncol=3,
                              nrow=2)
```

    ## 2024-05-07 17:14:30,009 - pegasusio.readwrite - INFO - zarr file '/projects/home/tlchan/projects/ascites/data_cite_objects/cd8.zarr.zip' is loaded.
    ## 2024-05-07 17:14:30,009 - pegasusio.readwrite - INFO - Function 'read_input' finished in 4.64s.

<img src="supp_cd8_figure_files/figure-gfm/fig_1A-1.png" width="1440" />

## Figure 1B

``` r
cd8_cluster_palette <- list("1" = "#FF0029",
                            "2" = "#377EB8",
                            "3" = "#66A61E",
                            "4" = "#984EA3",
                            "5" = "#00D2D5",
                            "6" = "#FF7F00",
                            "7" = "#AF8D00",
                            "8" = "#7F80CD",
                            "9" = "#B3E900",
                            "10" = "#C42E60",
                            "11" = "#A65628",
                            "12" = "#F781BF",
                            "13" = "#8DD3C7",
                            "14" = "#BEBADA",
                            "15" = "#FB8072"
)

cd8_data <- read.csv('/projects/home/tlchan/projects/ascites/figure_panels/abundance_data/cd8_patient_counts.csv') %>%
    mutate(Cluster = factor(Cluster))

ggplot(cd8_data, aes(x = Patient, y = Count, fill = Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_fill_manual(values = cd8_cluster_palette)
```

![](/tmp/supp_cd8_figure-4.rmd/supp_cd8_figure_files/figure-gfm/fig_1B-3.png)<!-- -->
