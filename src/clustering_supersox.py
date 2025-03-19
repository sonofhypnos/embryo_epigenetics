import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise
import matplotlib.pyplot as plt
import glob
import os

PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
SUPERSOX_RRBS_DATA_DIR = f"{PROJECT_DIR}data/supersox/rrbs/"
SUPERSOX_RNA_DATA = f"{PROJECT_DIR}data/supersox/rna/GSE247049_TPM1_hiPSCs_SOX2_vs_SOX2AV_vs_SOX2-17_hESCs_hFibs.csv"


def cluster_data(data, datatype="rrbs"):
    figname = f"{datatype}_cluster.png"
    # NOTE: I also tried:
    # data = data.dropna()
    # dist_matrix = pdist(data.T, metric="euclidean")
    dist_redundant = pairwise.nan_euclidean_distances(data.T)
    dist_matrix = squareform(dist_redundant)
    linkage_matrix = linkage(dist_matrix, method="average", metric="euclidean")

    # 4. Create visualization
    plt.figure(figsize=(12, 12))

    # Create dendrogram
    dendrogram(
        linkage_matrix,
        labels=data.columns,
        leaf_rotation=90,
        leaf_font_size=8,
    )

    plt.title(f"Hierarchical Clustering of {datatype} Data")
    plt.xlabel("Samples")
    plt.ylabel("Distance")
    plt.savefig(figname)
    plt.close()


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

    # extract names for columns
    basenames = [os.path.basename(f) for f in files]
    names = [extract_name(b) for b in basenames]

    # concatenate percent
    # TODO: what if we use the number of samples for the weights?
    full_dfs = [parse_data(f) for f in files]
    methylation_percent_series = [df["percent_methylated"] for df in full_dfs]
    df = pd.concat(methylation_percent_series, axis=1)

    df.columns = names

    cluster_data(df, datatype="rrbs")


def cluster_rna():
    rna = pd.read_csv(SUPERSOX_RNA_DATA)
    rna2 = rna.drop(["Name", "Gene ID"], axis=1)
    filtered_rna = rna2[
        (rna2 >= 1).all(axis=1)
    ]  # We tried .any, which also didn't result in in the original graph
    cluster_data(filtered_rna, datatype="rna")


cluster_rrbs()
cluster_rna()
