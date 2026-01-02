# spatialzones/plot.py
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

# ============================================================
#  Spatial and expression visualization utilities
# ============================================================

def plot_region_gene_spatial(temp, gene: str):
    """
    Plot spatial distribution of regions and gene expression.
    Uses region colors stored in temp.uns['region_colors'].

    Parameters
    ----------
    temp : AnnData
        AnnData object containing spatial data with `region` annotations.
    gene : str
        Gene name to visualize expression for.
    """
    # verify region colors exist
    if 'region_colors' not in temp.uns:
        raise ValueError("Missing temp.uns['region_colors']. Run assign_tumor_regions first.")
    if 'region' not in temp.obs:
        raise ValueError("Missing temp.obs['region']. Run assign_tumor_regions first.")

    sq.pl.spatial_scatter(
        temp,
        library_id="spatial",
        color=["region", gene],
        shape=None,
        size=0.5,
        img=False,
        ncols=2
    )


# ============================================================
#  Regional composition
# ============================================================

def plot_region_composition(temp):
    """
    Compute and plot the composition of cell types within each region.
    Uses stored region colors from temp.uns['region_colors'].
    """
    if 'region' not in temp.obs or 'temp_annotations' not in temp.obs:
        raise ValueError("Missing required obs columns: ['region', 'temp_annotations'].")

    region_order = list(temp.obs['region'].cat.categories)
    region_palette = dict(zip(region_order, temp.uns['region_colors']))

    comp = []
    for r in region_order:
        idx = np.where(temp.obs['region'] == r)[0]
        if len(idx) == 0:
            continue
        types = temp.obs.iloc[idx]['temp_annotations']
        cnt = types.value_counts()
        s = cnt / cnt.sum()
        comp.append(s.rename(r))

    if not comp:
        raise ValueError("No valid regions found in 'temp.obs['region']'.")

    comp_df = pd.concat(comp, axis=1).fillna(0)
    comp_df['Cell_Type'] = comp_df.index.astype(str)
    plot_df = comp_df.melt(id_vars='Cell_Type', var_name='Region', value_name='Fraction')

    plt.figure(figsize=(12, 6))
    sns.barplot(
        data=plot_df,
        x='Cell_Type',
        y='Fraction',
        hue='Region',
        hue_order=region_order,
        palette=region_palette
    )
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Fraction within region")
    plt.xlabel("Cell type")
    plt.title("Regional composition of cell types")
    plt.tight_layout()
    plt.show()


# ============================================================
#  Violin plot for gene by region
# ============================================================

def plot_gene_violin_by_region(temp, gene: str):
    """
    Violin plot of a gene's expression across regions using stored colors.
    """
    if 'region_colors' not in temp.uns:
        raise ValueError("Missing temp.uns['region_colors']. Run assign_tumor_regions first.")

    sc.pl.violin(
        temp,
        keys=gene,
        groupby='region',
        rotation=90,
        stripplot=False
    )


# ============================================================
#  Dotplot for gene by region
# ============================================================

def plot_gene_dotplot_by_region(temp, gene: str):
    """
    Dotplot of a gene's expression across regions.
    """
    sc.pl.dotplot(
        temp,
        var_names=[gene],
        groupby='region',
        cmap='viridis',
        standard_scale='var'
    )


# ============================================================
#  Graph vs Euclidean distance
# ============================================================

def plot_graph_vs_euclidean(temp):
    """
    Scatter plot comparing graph and Euclidean distance to tumor cells, colored by region.
    Adds Spearman correlation.
    """
    if not all(x in temp.obs.columns for x in ['graph_dist_to_tumor', 'euclid_dist_to_tumor']):
        raise ValueError("AnnData missing distance columns. Run assign_tumor_regions first.")

    region_order = list(temp.obs['region'].cat.categories)
    region_palette = dict(zip(region_order, temp.uns['region_colors']))

    rho, p = spearmanr(temp.obs['graph_dist_to_tumor'], temp.obs['euclid_dist_to_tumor'])

    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=temp.obs,
        x='graph_dist_to_tumor',
        y='euclid_dist_to_tumor',
        hue='region',
        hue_order=region_order,
        palette=region_palette,
        alpha=0.7,
        s=20
    )
    sns.regplot(
        data=temp.obs,
        x='graph_dist_to_tumor',
        y='euclid_dist_to_tumor',
        scatter=False,
        color='black',
        line_kws={'lw': 1, 'ls': '--', 'alpha': 0.7}
    )
    plt.title("Graph vs Euclidean distance to tumor")
    plt.text(
        0.05, 0.95,
        f"Spearman r = {rho:.2f}\np = {p:.1e}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
    )
    plt.tight_layout()
    plt.show()


# ============================================================
#  Real vs Expected distribution
# ============================================================

def plot_expected_vs_observed_regions(temp, remove_celltypes=None):
    """
    Side-by-side stacked barplots comparing:
      - Observed region composition per cell type
      - Expected (global) region composition

    Parameters
    ----------
    temp : AnnData
        Must contain temp.obs['region'], temp.obs['temp_annotations'],
        and temp.uns['region_colors']
    remove_celltypes : list of str or None
        Cell types to exclude (e.g. tumor cells)
    """
    if 'region' not in temp.obs:
        raise ValueError("Missing temp.obs['region']. Run assign_tumor_regions first.")
    if 'temp_annotations' not in temp.obs:
        raise ValueError("Missing temp.obs['temp_annotations'].")
    if 'region_colors' not in temp.uns:
        raise ValueError("Missing temp.uns['region_colors']. Run assign_tumor_regions first.")

    if remove_celltypes is None:
        remove_celltypes = []

    # Get region order and colors from the AnnData object
    region_order = list(temp.obs['region'].cat.categories)
    region_colors = dict(zip(region_order, temp.uns['region_colors']))

    # -----------------------
    # OBSERVED (per cell type)
    # -----------------------
    obs_df = temp.obs.copy()
    if remove_celltypes:
        obs_df = obs_df[~obs_df['temp_annotations'].isin(remove_celltypes)]

    obs_counts = (
        obs_df
        .groupby(['temp_annotations', 'region'])
        .size()
        .reset_index(name='count')
    )

    obs_counts['Proportion'] = (
        obs_counts
        .groupby('temp_annotations')['count']
        .transform(lambda x: x / x.sum())
    )

    obs_pivot = (
        obs_counts
        .pivot(index='temp_annotations', columns='region', values='Proportion')
        .fillna(0)
        [region_order]
    )

    # Order cell types by first region → last region gradient
    obs_pivot = obs_pivot.sort_values(
        by=[region_order[0], region_order[-1]], ascending=[False, True]
    )

    # -----------------------
    # EXPECTED (global)
    # -----------------------
    exp_counts = temp.obs['region'].value_counts()
    exp_prop = exp_counts / exp_counts.sum()

    exp_df = pd.DataFrame(
        {r: [exp_prop.get(r, 0)] for r in region_order},
        index=['Expected']
    )

    # -----------------------
    # PLOTTING
    # -----------------------
    fig, axes = plt.subplots(
        ncols=2,
        figsize=(18, 6),
        gridspec_kw={'width_ratios': [5, 0.6]}
    )

    # Observed plot
    obs_pivot.plot(
        kind='bar',
        stacked=True,
        color=[region_colors[r] for r in obs_pivot.columns],
        ax=axes[0],
        edgecolor='none'
    )

    axes[0].set_title("Observed regional composition per cell type", fontsize=13)
    axes[0].set_ylabel("Proportion")
    axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("")

    # --- FIXED TICKS ---
    axes[0].tick_params(axis='x', labelsize=10)
    axes[0].set_xticklabels(
        axes[0].get_xticklabels(),
        rotation=90,
        ha='right'
    )

    # Expected plot
    exp_df.plot(
        kind='bar',
        stacked=True,
        color=[region_colors[r] for r in exp_df.columns],
        ax=axes[1],
        edgecolor='none',
        legend=False
    )

    axes[1].set_title("Expected (global)", fontsize=13)
    axes[1].set_ylim(0, 1)
    axes[1].set_xlabel("")
    axes[1].set_xticklabels(['All cells'])

    # Legend
    axes[0].legend(
        title='Region',
        bbox_to_anchor=(1.02, 1),
        loc='upper left'
    )

    sns.despine(left=True, bottom=True)
    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout()
    plt.show()
