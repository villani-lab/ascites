library(tidyverse)
library(ggpubr)

plot_sf_boxplot <- function(analytes) {
    tissue_palette <- list("ascites" = "#00BFC4",
                           "plasma" = "#F8766D")

    sf_names <- read.csv("/projects/home/tlchan/projects/ascites/second_data_freeze/data/secreted_factors/sf_common_names.csv")
    paired_list <- list('ASC_41', 'ASC_43', 'ASC_45', 'ASC_46', 'ASC_48', 'ASC_52', 'ASC_57', 'ASC_61', 'ASC_62', 'ASC_65', 'ASC_66', 'ASC_67')

    # Load and prepare data
    sf_data <- read.csv("/projects/home/tlchan/projects/ascites/second_data_freeze/data/secreted_factors/secreted_factors_updated.csv", row.names = 1) %>%
        filter(diluted == "No") %>%
        filter(type != "pleural") %>%
        filter(panel == "96-cytokine") %>%
        filter(patient_id %in% paired_list) %>%
        mutate(log_concentration = log(concentration)) %>%
        merge(sf_names, all.x = TRUE) %>%
        mutate(analyte = ifelse(!is.na(common_name), common_name, analyte))

    # Select analytes
    sf_data <- sf_data %>%
        filter(analyte %in% analytes) %>%
        mutate(analyte = factor(analyte, levels = analytes))

    ggplot(sf_data, aes(x = type, y = log_concentration, fill = type)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.75) +
        geom_line(aes(group = patient_id)) +
        geom_point(pch = 20, size = 2) +
        stat_compare_means(paired = TRUE, label.x.npc = "center", aes(label = paste0("p = ", after_stat(p.format)))) +
        xlab("Type") +
        ylab("log(Concentration)") +
        facet_wrap(~analyte, scales = "free_y", nrow = 2) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
        scale_fill_manual(values = tissue_palette)
}