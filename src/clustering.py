import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
from sklearn.metrics import silhouette_score, pairwise
from scipy.cluster.hierarchy import fcluster

PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
WGBS_DATA_DIR = f"{PROJECT_DIR}data/"
RRBS_DATA_DIR = f"{PROJECT_DIR}rrbs_data/"
PBAT_DATA_DIR = f"{PROJECT_DIR}data/PBAT2017/"
SUPERSOX_RRBS_DATA_DIR = f"{PROJECT_DIR}data/supersox/rrbs/"
SUPERSOX_RNA_DATA = f"{PROJECT_DIR}data/supersox/rna/GSE247049_TPM1_hiPSCs_SOX2_vs_SOX2AV_vs_SOX2-17_hESCs_hFibs.csv"
CACHE_DIR = f"{PROJECT_DIR}tmp/"
RESULTS_DIR = f"{PROJECT_DIR}results/clustering/"


def create_if_not_exists(path):
    path = os.path.dirname(path)
    if not os.path.exists(path):
        os.mkdir(path)


create_if_not_exists(RESULTS_DIR)


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
    print(figure_filename)
    plt.close()

    return linkage_matrix


def analyze_cluster_separation(linkage_matrix, sample_groups):
    """
    Analyze the separation of sample groups in the clustering.

    Parameters:
    linkage_matrix: scipy linkage matrix from hierarchical clustering
    sample_groups: dictionary mapping sample names to their groups

    Returns:
    float: silhouette score measuring cluster quality
    """

    cluster_labels = fcluster(linkage_matrix, t=2, criterion="maxclust")

    # Calculate silhouette score
    group_labels = [sample_groups[sample] for sample in sample_names]
    score = silhouette_score(dist_matrix, group_labels)

    return score


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
        filtering = True
        if metric == "nan_euclidean":
            filtering = False
        linkage_matrix = cluster_data(df, metric=metric, filtering=filtering)


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
# cluster_rrbs_data(filtered_rna, metric="nan_euclidean", figure_dir=figure_dir)
