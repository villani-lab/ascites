library(ggpubr)
library(rstatix)
library(tidyverse)

plot_sf_boxplot <- function(analytes, nrow) {
    tissue_palette <- list("ascites" = "#2278B5",
                           "plasma" = "#D62A28")

    sf_names <- read.csv("/projects/home/tlchan/util/secreted_factors/sf_common_names.csv")
    paired_list <- list('ASC_41', 'ASC_43', 'ASC_45', 'ASC_46', 'ASC_48', 'ASC_52', 'ASC_57', 'ASC_61', 'ASC_62', 'ASC_65', 'ASC_66', 'ASC_67')

    # Load and prepare data
    sf_data <- read.csv("/projects/home/tlchan/projects/ascites/results/secreted_factors/integrated_data/secreted_factors_updated.csv", row.names = 1) %>%
        filter(diluted == "No") %>%
        filter(type != "pleural") %>%
        filter(panel == "96-cytokine") %>%
        filter(patient_id %in% paired_list) %>%
        mutate(log_concentration = log(concentration)) %>%
        merge(sf_names, all.x = TRUE) %>%
        mutate(analyte = common_name) %>%
        mutate(type = factor(type, levels = c("plasma", "ascites")))

    # Select analytes
    sf_data <- sf_data %>%
        filter(analyte %in% analytes) %>%
        mutate(analyte = factor(analyte, levels = analytes))

    stats <- sf_data %>%
        group_by(analyte) %>%
        t_test(log_concentration ~ type, paired = TRUE) %>%
        adjust_pvalue(method = "fdr") %>%
        add_significance() %>%
        add_xy_position(x = "type") %>%
        mutate(text_color = ifelse(p.adj < 0.05, "black", "gray"))

    stats$clean_pvals <- sapply(stats$p, function(x) {
        if (x > 0.01) {
            return(as.character(round(x, 2)))
        } else if (x < 0.01 & x > 0.001) {
            return(as.character(round(x, 3)))
        } else {
            formatC(x, format = "e", digits = 0)
        }
    })

    ggplot(sf_data, aes(x = type, y = log_concentration)) +
        geom_boxplot(outlier.shape = NA, aes(fill = type)) +
        geom_line(aes(group = patient_id)) +
        geom_point(pch = 20, size = 2) +
        labs(fill = "Tissue type") +
        xlab("") +
        ylab("log(Concentration)") +
        facet_wrap(~analyte, scales = "free_y", nrow = nrow) +
        theme_classic(base_size = 20) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_fill_manual(values = tissue_palette) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
        stat_pvalue_manual(data = stats, label = "p = {clean_pvals}", size = 6, color = "text_color") +
        scale_color_manual(values = c("black" = "#000000", "gray" = "#808080"), guide = "none")
}