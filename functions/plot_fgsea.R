library(fgsea)

plot_fgsea <- function(fgsea_res, ranks, genes, lin, var, gs) {
    nes <- round(fgsea_res$NES[fgsea_res$pathway == gs], 3)
    pval <- round(fgsea_res$pval[fgsea_res$pathway == gs], 3)
    n_genes <- fgsea_res$size[fgsea_res$pathway == gs]

    rnk <- rank(-ranks)
    ord <- order(rnk)

    statsAdj <- ranks[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathway <- unname(as.vector(na.omit(match(genes, names(statsAdj)))))
    pathway <- sort(pathway)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))

    diff <- (max(tops) - min(bottoms)) / 8

    ggplot(toPlot, aes(x = x, y = y)) +
        geom_line(color = "blue") +
        geom_hline(yintercept = 0, colour = "black") +
        geom_segment(data = data.frame(x = pathway),
                     mapping = aes(x = x, y = -0.15,
                                   xend = x, yend = -0.25),
                     linewidth = 0.4) +
        scale_y_continuous(expand = c(0.05, 0.05)) +
        xlab("Rank") +
        ylab("Enrichment score") +
        geom_text(aes(label = "")) +
        annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 1000, y = 0.9) +
        annotate("text", label = glue("p-value : {pval}"), x = length(ranks) - 1000, y = 0.8) +
        annotate("text", label = glue("# genes : {n_genes}"), x = length(ranks) - 1000, y = 0.7) +
        ggtitle(glue("{lin}, {var}, {gs} signature")) +
        theme_classic(base_size = 12)
}