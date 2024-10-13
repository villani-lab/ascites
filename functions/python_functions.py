import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def plot_umap(lin_data, palette, color="Cluster", width=6, height=6, legend_loc="on data", size=None):
    if size == None:
        size = 120000 / lin_data.shape[0]

    fig, ax = plt.subplots(1)
    sc.pl.umap(adata=lin_data.to_anndata(),
               color=color,
               use_raw=True,
               palette=palette,
               size=size,
               frameon=False,
               legend_loc=legend_loc,
               legend_fontoutline=5,
               title="",
               show=False,
               ax=ax)

    fig = plt.gcf()
    fig.set_size_inches(width, height)
    fig.tight_layout()
    ax.set_rasterization_zorder(2)

    return fig

def plot_feature(lin_data, genes, ncol, nrow):
    umap_coords = pd.DataFrame(lin_data.obsm['X_umap'], columns=['x', 'y'])
    x = umap_coords['x']
    y = umap_coords['y']

    # Set up figure structure
    fig_size = (5 * ncol, 4 * nrow)
    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=fig_size,
                             sharex=True, sharey=True)

    try:
        ax = axes.ravel()
    except:
        ax = [axes]

    # Plot for each gene
    for num, gene in enumerate(genes):
        # Get counts for each gene
        norm_counts = lin_data[:, gene].copy().get_matrix('X').todense().transpose()
        norm_counts = np.squeeze(np.asarray(norm_counts))

        if norm_counts.size == 0:  # ie. All empty/zero in sparce matrix
            norm_counts = np.zeros(norm_counts.shape[1])

        # Get n cells and perc cells
        ncells = (norm_counts != 0).sum()
        pcells = ncells / norm_counts.shape[0]

        if gene.startswith('cite_'):
            cmap = 'PuBu'
        else:
            cmap = 'YlOrRd'

        # Create heatmap
        hb = ax[num].hexbin(x=x,
                            y=y,
                            C=norm_counts,
                            cmap=cmap,
                            gridsize=150,
                            edgecolors="none")

        # Add percent expression
        ax[num].annotate(f'{ncells:,} ({pcells:.1%}) cells',
                         xy=(0.01, 0), xycoords='axes fraction',
                         fontsize=10,
                         horizontalalignment='left',
                         verticalalignment='bottom')

        # Add axes and colorbar information
        cb = fig.colorbar(hb, ax=ax[num], shrink=.75, aspect=10)
        cb.ax.set_title('logCPM', loc='left', fontsize=14)
        ax[num].set_title(gene, fontsize=18)
        ax[num].tick_params(left=False, labelleft=False,
                            bottom=False, labelbottom=False)

        if (num + ncol) % ncol == 0:
            # Start of row
            ax[num].set_ylabel('UMAP2', fontsize=18)
        if nrow == 1:
            # Only one row
            ax[num].set_xlabel('UMAP1', fontsize=18)
        elif num > (len(genes) - ncol - 1):
            # Last row if more than one row
            ax[num].set_xlabel('UMAP1', fontsize=18)

        ax[num].set_rasterization_zorder(2)

    for i in range(len(ax) - (len(ax) - len(genes)), len(ax)):
        ax[i].set_axis_off()
    fig.tight_layout()

    return fig


def plot_feature_by_tissue_type(lin_data, genes):
    # Get UMAP coordinates for base
    blood_data = lin_data[lin_data.obs['tissue_type'] == 'blood'].copy()
    blood_coords = pd.DataFrame(blood_data.obsm['X_umap'], columns=['x', 'y'])
    blood_x = blood_coords['x']
    blood_y = blood_coords['y']

    ascites_data = lin_data[lin_data.obs['tissue_type'] == 'ascites'].copy()
    ascites_coords = pd.DataFrame(ascites_data.obsm['X_umap'], columns=['x', 'y'])
    ascites_x = ascites_coords['x']
    ascites_y = ascites_coords['y']

    x = [blood_x, ascites_x]
    y = [blood_y, ascites_y]

    # set up figure structure
    ncol = 2
    nrow = len(genes)
    fig_size = (5 * ncol, 4 * nrow)
    fig, axes = plt.subplots(nrows=nrow, ncols=ncol, figsize=fig_size, sharex=True, sharey=True)
    ax = axes.ravel()
    title = ['Blood', 'Ascites']

    # Plot for each gene
    for num, gene in enumerate(genes):
        # Get gene counts
        blood_cts = blood_data[:, gene].copy().get_matrix('X').todense().transpose()
        blood_cts = np.squeeze(np.asarray(blood_cts))
        if blood_cts.size == 0:  # ie. All empty/zero in sparce matrix
            blood_cts = np.zeros(blood_cts.shape[1])
        blood_n = (blood_cts != 0).sum()
        blood_perc = blood_n / blood_cts.shape[0]

        ascites_cts = ascites_data[:, gene].copy().get_matrix('X').todense().transpose()
        ascites_cts = np.squeeze(np.asarray(ascites_cts))
        if ascites_cts.size == 0:  # ie. All empty/zero in sparce matrix
            ascites_cts = np.zeros(ascites_cts.shape[1])
        ascites_n = (ascites_cts != 0).sum()
        ascites_perc = (ascites_cts != 0).sum() / ascites_cts.shape[0]

        norm_counts = [blood_cts, ascites_cts]
        cb_max = max(max(blood_cts), max(ascites_cts))  # Makes sure colorbars are the same

        if gene.startswith('cite_'):
            cmap = 'PuBu'
        else:
            cmap = 'YlOrRd'

        # Create heatmap
        for i in range(2):
            hb = ax[num * 2 + i].hexbin(x=x[i],
                                        y=y[i],
                                        C=norm_counts[i],
                                        cmap=cmap,
                                        gridsize=150,
                                        vmin=0,
                                        vmax=cb_max,
                                        edgecolors="none")
            cb = fig.colorbar(hb, ax=ax[num * 2 + i], shrink=.75, aspect=10)
            cb.ax.set_title('logCPM', loc='left', fontsize=14)
            if i == 0:
                # ie. Blood
                ax[num * 2 + i].set_ylabel('UMAP2', fontsize=18)
                ax[num * 2 + i].annotate(f'{blood_n:,} ({blood_perc:.1%}) cells',
                                         xy=(0.01, 0), xycoords='axes fraction',
                                         fontsize=10,
                                         horizontalalignment='left',
                                         verticalalignment='bottom')
            if num + 1 == len(genes):
                # ie. Blood and only one row
                ax[num * 2 + i].set_xlabel('UMAP1', fontsize=18)
            if i == 1:
                # ie. Ascites
                ax[num * 2 + i].annotate(f'{ascites_n:,} ({ascites_perc:.1%}) cells',
                                         xy=(0.01, 0), xycoords='axes fraction',
                                         fontsize=10,
                                         horizontalalignment='left',
                                         verticalalignment='bottom')
            ax[num * 2 + i].set_title(f'{title[i]}: {gene}', fontsize=18)
            ax[num * 2 + i].tick_params(left=False, labelleft=False,
                                        bottom=False, labelbottom=False)
            ax[num * 2 + i].set_rasterization_zorder(2)

    fig.tight_layout()

    return fig


def plot_rss(rss, cell_type, top_n=5, max_n=None, ax=None):
    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(4, 4))
    if max_n is None:
        max_n = rss.shape[1]
    data = rss.T[cell_type].sort_values(ascending=False)[0:max_n]
    ax.plot(np.arange(len(data)), data, ".")
    ax.set_ylim([math.floor(data.min() * 100.0) / 100.0, math.ceil(data.max() * 100.0) / 100.0])
    ax.set_ylabel("RSS")
    ax.set_xlabel("Regulon")
    ax.set_title(cell_type)
    ax.set_xticklabels([])

    font = {
        "color": "red",
        "weight": "normal",
        "size": 14,
    }

    for idx, (regulon_name, rss_val) in enumerate(
            zip(data[0:top_n].index, data[0:top_n].values)
    ):
        ax.plot([idx, idx], [rss_val, rss_val], "r.")
        ax.text(
            idx + (max_n / 25),
            rss_val,
            regulon_name,
            fontdict=font,
            horizontalalignment="left",
            verticalalignment="center",
        )
