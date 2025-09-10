import numpy as np
import pandas as pd
import scanpy as sc
from anticor_features.anticor_features import get_anti_cor_genes
# from correlation import *
from sklearn import preprocessing
from sklearn.metrics import jaccard_score
from scipy import stats

def analyze_gene_correlations(expression_matrix, n_bins):
    """
    Analyze gene correlations in single-cell RNA expression data mimicing DubStepR.
    
    Parameters:
    expression_matrix (pd.DataFrame): Genes as rows, cells as columns
    
    Returns:
    list: Gene names that meet the correlation range and z-score criteria
    """
    
    def calculate_correlation_range(corr_matrix):
        """Calculate correlation range for each gene using 3rd highest to lowest difference"""
        ranges = {}
        for gene in corr_matrix.index:
            # Get correlations for this gene with all others
            correlations = corr_matrix[gene].drop(gene)  # Remove self-correlation
            # Sort correlations in descending order
            sorted_corrs = sorted(correlations, reverse=True)
            # Calculate range between 3rd highest and minimum
            if len(sorted_corrs) >= 3:
                corr_range = sorted_corrs[2] - sorted_corrs[-1]
                ranges[gene] = corr_range
        return pd.Series(ranges)
    
    def bin_by_expression(expression_matrix, n_bins=10):
        """Bin genes based on total expression"""
        total_expression = expression_matrix.sum(axis=1)
        log_expression = np.log1p(total_expression)
        # Create bins based on log expression
        bins = pd.qcut(log_expression, q=n_bins, labels=False)
        return bins
    
    # 1. Calculate gene-gene correlation matrix
    # Nomalize data
    data = pd.DataFrame.to_numpy(expression_matrix)
    data = preprocessing.normalize(data, norm='l1')
    normalized_df = pd.DataFrame(data)
    corr_matrix = pd.DataFrame(
        np.corrcoef(normalized_df.T),
        index=expression_matrix.columns,
        columns=expression_matrix.columns
    )
    
    # 2. Calculate correlation range for each gene
    corr_ranges = calculate_correlation_range(corr_matrix)
    
    # 3. Bin genes by expression level
    gene_bins = bin_by_expression(expression_matrix.T, n_bins)
    
    # 4. Calculate z-scores within each bin and filter genes
    selected_genes = []
    for bin_num in range(gene_bins.max() + 1):
        # Get genes in this bin
        bin_genes = gene_bins[gene_bins == bin_num].index
        if len(bin_genes) > 1:  # Need at least 2 genes to calculate z-score
            # Get correlation ranges for genes in this bin
            bin_ranges = corr_ranges[bin_genes]
            # Calculate z-scores
            z_scores = stats.zscore(bin_ranges)
            # Select genes with z-score > 0.7
            selected_genes.extend(
                bin_genes[z_scores > 0.7].tolist()
            )
    
    return selected_genes

def plot_cell_gene_heatmap(expression_matrix, highlighted_genes, 
                           title='Cell-Gene Expression Heatmap\n(Highlighted Genes in Red)',
                          figsize=(12, 8), 
                          cmap="coolwarm"):
    """
    Plot cell-gene expression heatmap with highlighted gene columns.
    
    Parameters:
    expression_matrix (pd.DataFrame): Expression matrix (cells as rows, genes as columns)
    highlighted_genes (list): List of gene names to highlight
    figsize (tuple): Figure size
    cmap (str): Colormap for the heatmap
    
    Returns:
    matplotlib.figure.Figure: The generated figure
    """
    # Create a copy of the data for plotting
    plot_data = expression_matrix.copy()
    
    # Log transform the data for better visualization
    plot_data = np.log1p(plot_data)
    
    # Create the figure and axes
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the heatmap
    sns.heatmap(plot_data, 
                cmap=cmap,
                yticklabels=False,  # Hide cell labels
                xticklabels=True)   # Show gene labels
    
    # Highlight the selected genes
    for gene in highlighted_genes:
        if gene in plot_data.columns:
            # Get the column index for the gene
            col_idx = plot_data.columns.get_loc(gene)
            
            # Add colored rectangle to highlight the column
            ax.add_patch(plt.Rectangle(
                (col_idx, 0),  # Starting point (x, y)
                1,            # Width
                plot_data.shape[0],  # Height (number of cells)
                fill=False,
                edgecolor='yellow',
                linewidth=2,
                zorder=10
            ))
            
            # Color and style the gene label
            label = ax.get_xticklabels()[col_idx]
            label.set_color('green')
            label.set_weight('bold')
            label.set_rotation(45)
            label.set_ha('right')
    
    plt.title(title, pad=20)
    plt.xlabel('Genes')
    plt.ylabel('Cells')
    
    # Adjust layout to prevent label cutoff
    plt.tight_layout()
    plt.show()

# 二值化
def binarize(data):
    """
    将NumPy数组二值化。
    非零值转换为1，零值保持为0。
    """
    return (data != 0).astype(int)

# 对于二值化的数据
def compute_jaccard_index(matrix):
    n_cols = matrix.shape[1]
    jaccard_matrix = np.zeros((n_cols, n_cols))
    
    for i in range(n_cols):
        for j in range(i, n_cols):
            jaccard_matrix[i, j] = jaccard_score(matrix[:, i], matrix[:, j])
            jaccard_matrix[j, i] = jaccard_score(matrix[:, i], matrix[:, j])
    
    return jaccard_matrix