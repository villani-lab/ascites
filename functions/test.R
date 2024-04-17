library(ggplot2)
library(glue)
library(rstatix)
library(tidyverse)
library(ggrepel)
library(xlsx)
library(ggpubr)

tissue_palette <- list("ascites" = "#00BFC4",
                       "plasma" = "#F8766D")

data_dir <- "/projects/home/tlchan/projects/ascites/second_data_freeze/data/secreted_factors"
fig_dir <- "/projects/home/tlchan/projects/ascites/second_data_freeze/figures/secreted_factors"

# Load and prepare data
sf_data <- read.csv(glue("{data_dir}/secreted_factors_updated.csv"), row.names = 1)
sf_data <- sf_data[sf_data["diluted"] == "No",]
sf_data <- sf_data[sf_data["type"] != "pleural",]
sf_data$log_concentration <- log(sf_data$concentration)

sf_names <- read.csv("/projects/home/tlchan/projects/ascites/second_data_freeze/data/secreted_factors/sf_common_names.csv")
sf_data <- merge(sf_data, sf_names, all.x = TRUE)
sf_data$analyte <- ifelse(!(is.na(sf_data$common_name)), sf_data$common_name, sf_data$analyte)


analyte_list <- c("CXCL10", "CXCL9", "CXCL16", "SDF-1", "IL-15", "IL-12p40", "TNFα", "IFNγ")
paired_list <- list('ASC_41', 'ASC_43', 'ASC_45', 'ASC_46', 'ASC_48', 'ASC_52', 'ASC_57', 'ASC_61', 'ASC_62', 'ASC_65', 'ASC_66', 'ASC_67')
sf_data <- sf_data %>%
    filter(panel == "96-cytokine") %>%
    filter(patient_id %in% paired_list) %>%
    filter(analyte %in% analyte_list) %>%
    mutate(analyte = factor(analyte, levels = analyte_list))

ggplot(sf_data, aes(x = type, y = log_concentration, fill = type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_line(aes(group = patient_id)) +
    geom_point(pch = 20, size = 2) +
    stat_compare_means(paired = TRUE, label.x.npc = "center", aes(label = paste0("p = ", after_stat(p.format)))) +
    xlab("Type") +
    ylab("log(Concentration)") +
    ggtitle("Gene concentration by tissue type") +
    facet_wrap(~analyte, scales = "free_y", nrow = 2) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_fill_manual(values = tissue_palette)