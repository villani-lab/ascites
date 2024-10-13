library(ggpubr)
library(parameters)
library(tidyverse)

plot_cluster_abundance <- function(lin, cluster_order, remove_clusters, n_breaks = 5) {
    tissue_palette <- list("ascites" = "#2278B5",
                           "blood" = "#D62A28",
                           "other" = "#000000")

    paired_samples <- c("ASC_10", "ASC_25", "ASC_41", "ASC_45", "ASC_46", "ASC_48", "ASC_49", "ASC_52", "ASC_57", "ASC_61", "ASC_62", "ASC_65", "ASC_66", "ASC_67")

    # Load data
    abundance <- read.csv("/projects/home/tlchan/projects/ascites/second_data_freeze/data/metadata/ascites_abundance.csv")

    # Remove cancer cells and doublets
    abundance <- abundance %>%
        filter(lineage != "cancer") %>%
        filter(patient_id %in% paired_samples) %>%
        filter(patient_id != "ASC_48")

    # Get immune count
    abundance <- abundance %>%
        group_by(patient_id, tissue_type, lineage, cluster) %>%
        mutate(clust_count = sum(immune == "True")) %>%
        select(patient_id, tissue_type, lineage, cluster, clust_count) %>%
        distinct()

    # Account for patient/tissue pairs that contributed zero cells to clusters
    combos <- expand.grid(unique(abundance$patient_id), unique(abundance$tissue_type), unique(abundance$lineage), unique(abundance$cluster)) %>% `colnames<-`(c("patient_id", "tissue_type", "lineage", "cluster"))

    abundance <- combos %>%
        left_join(abundance, by = c("patient_id", "tissue_type", "lineage", "cluster")) %>%
        replace(is.na(.), 0)

    # Get statistics
    abundance <- abundance %>%
        group_by(patient_id, tissue_type) %>%
        mutate(total_count = sum(clust_count)) %>%
        mutate(clust_percentage = clust_count / total_count * 100) %>%
        mutate(log_clust_percentage = log1p(clust_percentage)) %>%
        mutate(tissue_type = factor(tissue_type, levels = c("ascites", "blood")))

    # Remove patients with fewer than 250 immune native fraction cells
    abundance <- abundance %>%
        group_by(patient_id) %>%
        mutate(min_total_count = min(total_count)) %>%
        filter(min_total_count > 250)

    # Subset abundance
    lin_abundance <- abundance %>% filter(lineage == lin)

    # Remove extra clusters
    lin_abundance <- lin_abundance %>%
        group_by(cluster) %>%
        filter(sum(clust_count) != 0) %>%
        droplevels()

    # Order clusters
    if (hasArg(cluster_order)) {
        lin_abundance <- lin_abundance %>% mutate(cluster = factor(cluster, levels = cluster_order))
    } else {
        lin_abundance <- lin_abundance %>% mutate(cluster = factor(cluster))
    }

    # Remove selected clusters
    if (hasArg(remove_clusters)) {
        lin_abundance <- lin_abundance %>% filter(!cluster %in% remove_clusters)
    }

    bp <- ggplot(lin_abundance, aes(x = clust_percentage + 1, y = cluster, fill = tissue_type)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(pch = 21, position = position_jitterdodge(), aes(fill = tissue_type), size = 2) +
        scale_x_log10() +
        coord_cartesian(clip = "off") +
        scale_y_discrete(limits = rev) +
        labs(fill = "Tissue type") +
        xlab("Percent native immune + 1") +
        ylab("") +
        theme_classic(base_size = 27) +
        theme(axis.text.y = element_blank(), axis.text = element_text(size = 20)) +
        scale_fill_manual(values = tissue_palette)

    pt_res <- lapply(unique(lin_abundance$cluster), function(clust) {
        pt_data <- lin_abundance %>% filter(cluster == clust)
        stats <- parameters(t.test(log_clust_percentage ~ tissue_type, pt_data, paired = TRUE))
        stats$cluster <- clust
        return(stats)
    }) %>%
        do.call(rbind, .) %>%
        mutate(padj = p.adjust(p, method = "fdr")) %>%
        mutate(color = case_when(padj < 0.1 & Difference > 0 ~ "ascites", padj < 0.1 & Difference < 0 ~ "blood", padj >= 0.1 ~ "other"))

    fp <- ggplot(pt_res, aes(x = Difference, y = factor(cluster), color = color)) +
        geom_point(size = 3) +
        geom_errorbarh(mapping = aes(xmin = CI_low, xmax = CI_high, height = 0)) +
        geom_vline(xintercept = 0) +
        scale_y_discrete(limits = rev) +
        guides(color = "none") +
        xlab("Log2FC") +
        ylab("") +
        scale_x_continuous(n.breaks = n_breaks) +
        theme_classic(base_size = 27) +
        theme(axis.text = element_text(size = 20)) +
        scale_color_manual(values = tissue_palette)

    ggarrange(fp, bp, ncol = 2, nrow = 1, widths = c(0.5, 1.0), common.legend = TRUE, legend = "bottom")
}