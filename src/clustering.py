import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
from sklearn.metrics import silhouette_score, pairwise
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform

PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
WGBS_DATA_DIR = f"{PROJECT_DIR}data/"
RRBS_DATA_DIR = f"{PROJECT_DIR}rrbs_data/"
PBAT_DATA_DIR = f"{PROJECT_DIR}data/PBAT2017/"
SUPERSOX_RRBS_DATA_DIR = f"{PROJECT_DIR}data/supersox/rrbs/"
SUPERSOX_RNA_DATA = f"{PROJECT_DIR}data/supersox/rna/GSE247049_TPM1_hiPSCs_SOX2_vs_SOX2AV_vs_SOX2-17_hESCs_hFibs.csv"
CACHE_DIR = f"{PROJECT_DIR}tmp/"
RESULTS_DIR = f"{PROJECT_DIR}results/clustering/"


def prepare_weighted_data(files):
    """
    Returns methylation values and their corresponding coverages
    """
    full_dfs = [parse_data(f) for f in files]
    methylation_series = [df["percent_methylated"] for df in full_dfs]
    coverage_series = [df["meth_count"] + df["unmeth_count"] for df in full_dfs]

    methylation_df = pd.concat(methylation_series, axis=1)
    coverage_df = pd.concat(coverage_series, axis=1)

    names = [extract_name(os.path.basename(f)) for f in files]
    methylation_df.columns = names
    coverage_df.columns = names

    return methylation_df, coverage_df


def beta_variance_weight(n, max_weight=1.0):
    """
    Convert coverage to weights based on beta distribution variance

    Parameters:
    methylated: number of methylated reads
    unmethylated: number of unmethylated reads
    max_weight: scaling factor for maximum weight
    """
    weight = n / (n + 1)  # approaches 1 as n increases
    return max_weight * weight


def weighted_distance(methylation_data, coverage_data, metric="euclidean"):
    """
    Compute pairwise distances with coverage weights
    """
    n_samples = methylation_data.shape[1]
    dist_matrix = np.zeros((n_samples, n_samples))

    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            # Get values and coverages for pair of samples
            values1 = methylation_data.iloc[:, i]
            values2 = methylation_data.iloc[:, j]
            weights1 = coverage_data.iloc[:, i]
            weights2 = coverage_data.iloc[:, j]

            # Use only positions where both samples have data
            valid = ~(values1.isna() | values2.isna())
            if not valid.any():
                dist_matrix[i, j] = dist_matrix[j, i] = np.nan
                continue

            v1 = values1[valid]
            v2 = values2[valid]
            w1 = beta_variance_weight(weights1[valid])
            w2 = beta_variance_weight(weights2[valid])

            # Combined weight for each position
            combined_weights = np.minimum(w1, w2)  # or could use mean/harmonic mean

            if metric == "euclidean":
                # Weighted euclidean distance
                diff_squared = (v1 - v2) ** 2
                dist = np.sqrt(np.average(diff_squared, weights=combined_weights))

            elif metric == "correlation":
                # Weighted correlation
                mean1 = np.average(v1, weights=combined_weights)
                mean2 = np.average(v2, weights=combined_weights)
                centered1 = v1 - mean1
                centered2 = v2 - mean2
                numerator = np.average(centered1 * centered2, weights=combined_weights)
                var1 = np.average(centered1**2, weights=combined_weights)
                var2 = np.average(centered2**2, weights=combined_weights)
                correlation = numerator / np.sqrt(var1 * var2)
                dist = 1 - correlation

            elif metric == "angle":
                # Weighted cosine distance
                numerator = np.average(v1 * v2, weights=combined_weights)
                norm1 = np.sqrt(np.average(v1**2, weights=combined_weights))
                norm2 = np.sqrt(np.average(v2**2, weights=combined_weights))
                cosine = numerator / (norm1 * norm2)
                dist = np.arccos(np.clip(cosine, -1.0, 1.0))

            dist_matrix[i, j] = dist_matrix[j, i] = dist

    return squareform(dist_matrix)


def create_if_not_exists(path):
    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.makedirs(path)


create_if_not_exists(RESULTS_DIR)


# Modified cluster_data function:
def cluster_data_weighted(
    methylation_data,
    coverage_data=None,
    var_threshold=0.1,
    metric="correlation",
    figure_dir=RESULTS_DIR,
):
    filtering = True
    if metric == "nan_euclidean":
        filtering = False
    filter_string = "filtered/" if filtering else "unfiltered/"
    figure_filename = f"{figure_dir}{filter_string}hierarchy_{metric}_weighted.png"

    if filtering:
        valid_rows = ~methylation_data.isna().any(axis=1)
        methylation_data = methylation_data[valid_rows]
        if coverage_data is not None:
            coverage_data = coverage_data[valid_rows]

    # # Filter for variable sites
    # site_vars = methylation_data.var(axis=1)
    # variable_sites = methylation_data[site_vars > var_threshold]
    # if coverage_data is not None:
    #     variable_coverage = coverage_data[site_vars > var_threshold]

    # Calculate distance matrix
    if coverage_data is not None:
        dist_matrix = weighted_distance(methylation_data, coverage_data, metric=metric)
    else:
        dist_matrix = pdist(methylation_data.T, metric=metric)
    # print(dist_matrix)
    # 3. Perform hierarchical clustering
    linkage_matrix = linkage(dist_matrix, method="euclidean")

    create_if_not_exists(figure_filename)
    # 4. Create visualization
    plt.figure(figsize=(12, 12))

    # Create dendrogram
    dendrogram(
        linkage_matrix,
        labels=methylation_data.columns,
        leaf_rotation=90,
        leaf_font_size=8,
    )

    plt.title(f"Hierarchical Clustering of RRBS Methylation Data ({metric})")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.savefig(figure_filename)
    plt.close()

    return linkage_matrix


def cluster_data(
    methylation_data,
    sample_metadata=None,
    var_threshold=0.1,
    metric="correlation",
    figure_dir=RESULTS_DIR,
):
    filtering = True
    if metric == "nan_euclidean":
        filtering = False
    filter_string = "filtered/" if filtering else "unfiltered/"
    figure_filename = f"{figure_dir}{filter_string}hierarchy_{metric}.png"
    """
    Perform hierarchical clustering on RRBS methylation data.

    Parameters:
    methylation_data (pd.DataFrame): DataFrame with CpG sites as rows, samples as columns,
                                   values are methylation beta values (0-1)
    sample_metadata (pd.DataFrame): Optional metadata about samples for annotation

    Returns:
    fig: matplotlib figure object
    linkage_matrix: scipy linkage matrix
    """
    # 1. Data preprocessing
    # Remove CpG sites with missing values

    print(metric)

    if filtering:
        methylation_data = methylation_data.dropna()

    # Optional: Filter for variable CpG sites
    site_vars = methylation_data.var(axis=1)
    variable_sites = methylation_data[site_vars > var_threshold]

    sites = variable_sites

    if metric == "angle":
        dist_matrix = pdist(variable_sites.T, metric="cosine")
        dist_matrix = np.arccos(
            1 - dist_matrix
        )  # Convert cosine distance to angular distance
    elif metric == "nan_euclidean":
        dist_redundant = pairwise.nan_euclidean_distances(variable_sites.T)
        dist_matrix = squareform(dist_redundant)
    else:
        # 2. Calculate distance matrix
        # Using correlation distance (1 - correlation)
        dist_matrix = pdist(methylation_data.T, metric=metric)

    # print(dist_matrix)
    # 3. Perform hierarchical clustering
    if metric == "euclidean":
        linkage_matrix = linkage(dist_matrix, method="single", metric="euclidean")
    else:
        linkage_matrix = linkage(dist_matrix, method="ward")

    create_if_not_exists(figure_filename)
    # 4. Create visualization
    plt.figure(figsize=(12, 12))

    # Create dendrogram
    dendrogram(
        linkage_matrix,
        labels=sites.columns,
        leaf_rotation=90,
        leaf_font_size=8,
    )

    plt.title("Hierarchical Clustering of RRBS Methylation Data")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.savefig(figure_filename)
    plt.close()

    return linkage_matrix


def parse_data(filename):
    methylation = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        comment="#",
    )

    methylation.columns = [
        "chrom",
        "start",
        "end",
        "percent_methylated",
        "meth_count",
        "unmeth_count",
        "DMR",
    ]
    methylation = methylation.set_index("start")
    return methylation


def extract_name(filename):
    # Remove the end part
    without_end = filename.split("_DMRs_bvals.bed")[0]
    # Remove the GSM number
    without_gsm = without_end.split("GSM")[1].split("_", 1)[1]
    return without_gsm


def cluster_rrbs():
    files = glob.glob(f"{SUPERSOX_RRBS_DATA_DIR}*.bed")
    basenames = [os.path.basename(f) for f in files]
    names = [extract_name(b) for b in basenames]

    full_dfs = [parse_data(f) for f in files]
    methylation_percent_series = [df["percent_methylated"] for df in full_dfs]
    df = pd.concat(methylation_percent_series, axis=1)
    df.columns = names

    for metric in ["angle", "correlation", "euclidean", "nan_euclidean"]:
        linkage_matrix = cluster_data(df, metric=metric)


def cluster_rrbs_weighted():
    files = glob.glob(f"{SUPERSOX_RRBS_DATA_DIR}*.bed")
    methylation_df, coverage_df = prepare_weighted_data(files)

    cluster_data_weighted(
        methylation_df,
        coverage_df,
        metric="euclidean",
        figure_dir=f"{RESULTS_DIR}weighted/",
    )


def cluster_rna():
    rna = pd.read_csv(SUPERSOX_RNA_DATA)

    rna2 = rna.drop(["Name", "Gene ID"], axis=1)
    filtered_rna = rna2[(rna2 >= 1).all(axis=1)]
    print(filtered_rna)
    figure_dir = f"{RESULTS_DIR}rna_filtered/"
    create_if_not_exists(figure_dir)
    cluster_data(filtered_rna, metric="euclidean", figure_dir=figure_dir)
    figure_dir = f"{RESULTS_DIR}rna_unfiltered/"
    create_if_not_exists(figure_dir)

    cluster_data(rna2, metric="euclidean", figure_dir=figure_dir)


# cluster_rrbs_weighted()
cluster_rrbs()
# cluster_rrbs_data(filtered_rna, metric="nan_euclidean", figure_dir=figure_dir)
