import subprocess
import os
from itertools import product
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import diskcache
from diskcache import Cache
from diskcache.core import ENOVAL
from functools import partial
import glob
from dataclasses import dataclass, field
from typing import Any, Dict, List, Iterator
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
WGBS_DATA_DIR = f"{PROJECT_DIR}data/"
RRBS_DATA_DIR = f"{PROJECT_DIR}rrbs_data/"
SUPERSOX_DATA_DIR = f"{PROJECT_DIR}data/supersox/rrbs/"
CACHE_DIR = f"{PROJECT_DIR}tmp/"
RESULTS_DIR = f"{PROJECT_DIR}results/"

cache = Cache(CACHE_DIR)

# TODO: Add other methylation types from RRBS data
# TODO: Modify the GC data by taking some running average and see at which point the correlations actually get to a higher point (although there will probably not be a clear transition)
# TODO: Refactor plotting functions (add display names by default, add plt.close etc.)


CELL_TYPES_RRBS = [
    "MII_Oocyte+RRBS",
    "MII_Oocyte+scRRBS",
    "Sperm+RRBS",
    "Sperm+scRRBS",
    "Zygote+RRBS",
    "2-cell+RRBS",
    "4-cell+RRBS",
    "8-cell+RRBS",
    "Morula+RRBS",
    "ICM+RRBS",
    "ICM+WGBS",
    "TE+RRBS",
    "Liver+WGBS",
    "2nd_PB+RRBS",
    "1st_PB+RRBS",
    "PN+scRRBS",
]

CELL_TYPES_WGBS = [
    "Oocyte",
    "Sperm",
    "Zygote",
    "2cell",
    "4cell",
    "8cell",
    "Morula",
    "ICM",
    "TE",
    "hESC",
]
CELL_TYPES_SUPERSOX = [
    "O_Fibs_L2",
    "Y_Fibs_L2",
    "Y_iPSC_SS",
    "Y_iPSC_WT",
    "Y_iPSC_AV",
    "O_iPSC_SS",
    "O_iPSC_WT",
    "O_iPSC_AV",
    # NOTE: ignoring L2 for now
]

methylation_types_wgbs = ["ACG.TCG", "GCA.GCC.GCT"]
methylation_types_rrbs = ["CpG"]

# TODO: handle scRRBS data once I know what that is
# TODO: handle files where we have RRBS and WGBS like ICM?
cell_type_to_seq_type_rrbs = {
    "MII_Oocyte+RRBS": "RRBS",
    "MII_Oocyte+scRRBS": "scRRBS",
    "Sperm+RRBS": "RRBS",
    "Sperm+scRRBS": "scRRBS",
    "Zygote+RRBS": "RRBS",
    "2-cell+RRBS": "RRBS",
    "4-cell+RRBS": "RRBS",
    "8-cell+RRBS": "RRBS",
    "Morula+RRBS": "RRBS",
    "ICM+RRBS": "RRBS",  # WGBS also possible
    "ICM+WGBS": "WGBS",  # WGBS also possible
    "TE+RRBS": "RRBS",
    "Liver+WGBS": "WGBS",
    "2nd_PB+RRBS": "RRBS",
    "1st_PB+RRBS": "RRBS",
    "PN+scRRBS": "scRRBS",
}

cell_type_to_seq_type_wgbs = {cell_type: "WGBS" for cell_type in CELL_TYPES_WGBS}
cell_type_to_seq_type_supersox = {
    cell_type: "RRBS" for cell_type in CELL_TYPES_SUPERSOX
}


wgbs_filenames: Dict[str, Dict[str, str]] = {
    methylation_type: {
        cell_type: f"{WGBS_DATA_DIR}{cell_type}_{methylation_type}.wig"
        for cell_type in CELL_TYPES_WGBS
    }
    for methylation_type in methylation_types_wgbs
}

rrbs_filenames: Dict[str, Dict[str, str]] = {
    methylation_type: {
        cell_type: f"{RRBS_DATA_DIR}{cell_type.split('+')[0]}_{methylation_type}_{cell_type_to_seq_type_rrbs[cell_type]}.wig"
        for cell_type in CELL_TYPES_RRBS
    }
    for methylation_type in methylation_types_rrbs
}

supersox_filenames: Dict[
    str, Dict[str, str]
] = {  # FIXME: actually merge supersox files
    "RRBS": {
        cell_type: glob.glob(f"{SUPERSOX_DATA_DIR}GS*{cell_type}*.wig")[0]
        for cell_type in CELL_TYPES_SUPERSOX
    }
}

supersox_sample_filenames: Dict[str, Dict[str, str]] = {
    "RRBS": {
        cell_type: glob.glob(f"{SUPERSOX_DATA_DIR}GS*{cell_type}*.wig")
        for cell_type in CELL_TYPES_SUPERSOX
    }
}

# NOTE: I don't have the storage for this one (nor the data bandwidth)! (50GB)
wgbs_sample_filenames: Dict[str, Dict[str, str]] = {
    methylation_type: {
        cell_type: glob.glob(f"{WGBS_DATA_DIR}GS*{cell_type}*{methylation_type}*.bw")
        for cell_type in CELL_TYPES_WGBS
    }
    for methylation_type in methylation_types_wgbs
}

rrbs_sample_filenames: Dict[str, Dict[str, str]] = {
    methylation_type: {
        cell_type: glob.glob(
            f"{RRBS_DATA_DIR}GS*{cell_type_to_seq_type_rrbs[cell_type]}_{cell_type.split('+')[0]}*_methylation_calling_{methylation_type}.bw"
        )
        for cell_type in CELL_TYPES_RRBS
    }
    for methylation_type in methylation_types_rrbs
}


# print(wgbs_filenames)
# print(rrbs_filenames)
# print(rrbs_sample_filenames)
# print(wgbs_sample_filenames)


@dataclass
class DatasetParams:
    methylation_types: List[str]
    filenames: Dict[str, Dict[str, str]]
    sample_filenames: Dict[str, Dict[str, str]]
    cell_types: List[str]
    cell_type_seq_type: Dict[str, str]
    data_path: str
    results_path: str
    name: str
    # Allow dictionary-style access as a fallback
    def __getitem__(self, item):
        return getattr(self, item)


@dataclass
class AllDatasetParams:
    RRBS: DatasetParams
    WGBS: DatasetParams
    supersox: DatasetParams

    # Allow dictionary-style access as a fallback
    def __getitem__(self, item: str) -> DatasetParams:
        return getattr(self, item)

    def __iter__(self) -> Iterator[DatasetParams]:
        for dataset in (self.RRBS, self.WGBS, self.supersox):
            yield dataset


dataset_params = AllDatasetParams(
    RRBS=DatasetParams(
        methylation_types=methylation_types_rrbs,
        filenames=rrbs_filenames,
        sample_filenames=rrbs_sample_filenames,
        cell_types=CELL_TYPES_RRBS,
        cell_type_seq_type=cell_type_to_seq_type_rrbs,
        data_path=RRBS_DATA_DIR,
        results_path=os.path.join(RESULTS_DIR, "rrbs/"),
        name="RRBS",
    ),
    WGBS=DatasetParams(
        methylation_types=methylation_types_wgbs,
        filenames=wgbs_filenames,
        sample_filenames=wgbs_sample_filenames,
        cell_types=CELL_TYPES_WGBS,
        cell_type_seq_type=cell_type_to_seq_type_wgbs,
        data_path=WGBS_DATA_DIR,
        results_path=os.path.join(RESULTS_DIR, "wgbs/"),
        name="WGBS",
    ),
    supersox=DatasetParams(
        methylation_types=["RRBS"],
        filenames=supersox_filenames,
        sample_filenames=supersox_sample_filenames,
        cell_types=CELL_TYPES_SUPERSOX,
        cell_type_seq_type=cell_type_to_seq_type_supersox,
        data_path=SUPERSOX_DATA_DIR,
        results_path=os.path.join(RESULTS_DIR, "supersox/"),
        name="supersox",
    ),
)


def create_if_not_exists(path):
    if not os.path.exists(path):
        os.mkdir(path)


for dataset in dataset_params:
    create_if_not_exists(dataset.results_path)
    create_if_not_exists(dataset.data_path)


@cache.memoize()
def run_conda_command(
    command,
    env_name="epi_env",
    conda_path="/home/tassilo/miniconda3",
    capture_output=True,
):
    """
    Execute a command in a specific conda environment.

    Args:
        command (str): The command to execute
        env_name (str): Name of the conda environment
        conda_path (str): Path to conda installation
        capture_output (bool): Whether to capture and return command output

    Returns:
        tuple: (return_code, stdout, stderr) if capture_output=True
        subprocess.CompletedProcess: if capture_output=False
    """
    # Setup the conda initialization script
    conda_setup = f"""
    __conda_setup="$('{conda_path}/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "{conda_path}/etc/profile.d/conda.sh" ]; then
            . "{conda_path}/etc/profile.d/conda.sh"
        else
            export PATH="{conda_path}/bin:$PATH"
        fi
    fi
    unset __conda_setup

    conda activate {env_name}
    {command}
    """

    try:
        result = subprocess.run(
            ["zsh", "-c", conda_setup],
            capture_output=capture_output,
            text=True,
            check=False,
        )
        return result

    except subprocess.SubprocessError as e:
        if capture_output:
            return 1, "", str(e)
        raise


# Calculate figure size based on number of cell types
def get_figure_size(n_rows, n_cols):
    # Allow roughly 1 inch per cell type, plus margins
    width = max(10, n_cols * 1.2)  # minimum width of 10 inches
    height = max(8, n_rows * 1.2)  # minimum height of 8 inches
    return (width, height)


def correlation_to_angular_distance(r):
    """
    Convert correlation coefficient to angular distance in radians.

    Parameters:
    - r : float
        Correlation coefficient, must be in the range [-1, 1].

    Returns:
    - float
        Angular distance in radians.
    """
    r = round(r, 4)
    if r < -1 + 1e-6 or r - 1e-6 > 1:
        raise ValueError(
            f"Correlation coefficient must be between -1 and 1, value is {r}"
        )

    return np.arccos(r) / np.pi


# @cache.memoize()
def correlation(f1, f2, capture_output=True):
    if f1 == f2:
        return 1
    try:
        # results = subprocess.run(
        #     ["zsh", "-c", f"{PROJECT_DIR}src/wig_correlation {f1} {f2}"],
        #     capture_output=capture_output,
        #     text=True,
        #     check=False,
        # )

        results = run_conda_command(f"wigCorrelate {f1} {f2}")
        # print(results)
        corr = float(results.stdout.split()[-1])
    except TypeError as e:
        print(results)
        print(e)
    except Exception as e:
        print(e)
        raise e

    # print(result)

    return corr


def filepath_to_index(filepath):
    return os.path.basename(filepath).split("_")[0]


def generate_sample_histograms(params):
    for methylation_type in params.methylation_types:
        for cell_type in params.cell_types:
            for sample_filename in params.sample_filenames[methylation_type][cell_type]:
                hist_file = os.path.join(
                    params.results_path,
                    f"histogram_{filepath_to_index(sample_filename)}_{methylation_type}_{cell_type}.txt",
                )
                # print(hist_file)
                # print("exiting..")
                # exit()
                run_conda_command(
                    f"wiggletools histogram {hist_file} 100 {sample_filename}"
                )


def sample_correlations(dataset="RRBS", methylation_type="CpG"):
    params = dataset_params[dataset]
    for cell_type in params["cell_types"]:
        sample_filenames = params["sample_filenames"][methylation_type][cell_type]

        if len(sample_filenames) < 2:
            print(
                f"Not enough sample filenames for {cell_type} {methylation_type} {dataset} ({len(sample_filenames)} files)"
            )
            continue
        # if (
        #     cell_type == "8cell"
        # ):  # NOTE: 8cell files have some problem with their encoding leading to some issue with wigCorrelate
        #     continue

        def get_sample_correlations(sample_filenames):
            # if cell_type.split(" ")[0] in ["8cell", "4cell", "Morula"]:
            #     key = run_conda_command.__cache_key__(
            #         f"wigCorrelate {' '.join(sample_filenames)}"
            #     )

            #     print(f"Key: {key}")
            #     popped = cache.pop(key)
            #     print(f"Popped: ({popped})")

            result = run_conda_command(f"wigCorrelate {' '.join(sample_filenames)}")
            correlations = result.stdout

            samples = [sample.split("\t") for sample in correlations.split("\n")][
                :-1
            ]  # removing empty end
            # print(f"cell_type:{cell_type}")
            # print(f"samples:{samples}")

            files = list(set([file for sample in samples for file in sample[:2]]))

            # print(f"files:{files}")
            length = len(files)
            corr = np.zeros((length, length))
            for sample in samples:
                i = files.index(sample[0])
                j = files.index(sample[1])
                val = float(sample[2])
                corr[i, j] = val
                corr[j, i] = val

            return corr

        def get_sample_correlations_pair(sample_filenames):
            length = len(sample_filenames)
            corr = np.zeros((length, length))
            for i, filename1 in enumerate(sample_filenames):
                for j, filename2 in enumerate(sample_filenames):
                    pearson = correlation(filename1, filename2)
                    print(pearson)
                    corr[i, j] = pearson
            return corr

        indexes = [filepath_to_index(name) for name in sample_filenames]
        # if cell_type.split(" ")[0] in ["8cell", "4cell", "Morula"]:

        #     corr = get_sample_correlations_pair(sample_filenames)
        #     print(f"Cell-type: {cell_type}")
        #     print(corr)
        #     pause = input("Continue?")
        # else:
        corr = get_sample_correlations(sample_filenames)
        length = corr.shape[0]

        figsize = get_figure_size(length, length)

        filename = f"corr_matrix_{methylation_type}_{dataset}_{cell_type}.png"
        print(f"filename: {filename}")
        plt.figure(figsize=figsize)  # Adjust size as needed
        sns.heatmap(
            corr,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            xticklabels=indexes,
            yticklabels=indexes,
        )
        plt.title(
            f"Correlation between {cell_type} {methylation_type} samples in the {dataset} dataset"
        )
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(
            os.path.join(
                params.results_path,
                filename,
            )
        )
        plt.close()


def dataset_correlations(dataset_params=dataset_params):
    for dataset in dataset_params:
        params = dataset_params[dataset]
        for methylation_type in params["methylation_types"]:
            cell_types = params["cell_types"]
            file_names = params["filenames"]

            length = len(cell_types)
            corr = np.zeros((length, length))
            ang = np.zeros((length, length))
            figsize = get_figure_size(length, length)

            for i, cell1 in enumerate(cell_types):
                for j, cell2 in enumerate(cell_types):
                    if methylation_type == "GCA.GCC.GCT" and (
                        cell2 == "ICM" or cell1 == "ICM"
                    ):
                        continue
                        # FIXME: ICM file accidentally wasn't transferred to my pc (probably low disk space)!
                    cell_file_1 = file_names[methylation_type][cell1]
                    cell_file_2 = file_names[methylation_type][cell2]

                    corr[i, j] = correlation(cell_file_1, cell_file_2)
                    # print(cell1, cell2, corr[i, j])

                    ang[i, j] = correlation_to_angular_distance(corr[i, j])
            # print(corr)

            plt.figure(figsize=figsize)
            sns.heatmap(
                corr,
                annot=True,
                fmt=".2f",
                cmap="coolwarm",
                xticklabels=cell_types,
                yticklabels=cell_types,
            )
            plt.title(
                f"Correlation Matrix for {methylation_type} for {dataset} dataset"
            )
            plt.savefig(
                os.path.join(
                    params.results_path,
                    f"correlation_matrix_{methylation_type}_{dataset}.png",
                )
            )
            plt.close()


# For cell types that exist in both datasets, create display names
def get_display_name(cell_type, dataset):
    seq_type = dataset_params[dataset]["cell_type_seq_type"][cell_type]
    return f"{cell_type.split('+')[0]}_{seq_type}"


def inter_dataset_correlations():
    dataset_name_1 = "WGBS"
    dataset_name_2 = "RRBS"

    cell_types1 = dataset_params[dataset_name_1]["cell_types"]
    cell_types2 = dataset_params[dataset_name_2]["cell_types"]
    file_names1 = dataset_params[dataset_name_1]["filenames"]
    file_names2 = dataset_params[dataset_name_2]["filenames"]
    methylation_type_1 = "ACG.TCG"
    methylation_type_2 = "CpG"

    length1 = len(cell_types1)
    length2 = len(cell_types2)
    corr = np.zeros((length1, length2))
    ang = np.zeros((length1, length2))

    figsize = get_figure_size(length1, length2)

    # Create display names for axis labels
    display_names1 = [get_display_name(ct, dataset_name_1) for ct in cell_types1]
    display_names2 = [get_display_name(ct, dataset_name_2) for ct in cell_types2]

    for i, cell1 in enumerate(cell_types1):
        for j, cell2 in enumerate(cell_types2):
            if methylation_type_1 == "GCA.GCC.GCT" and (
                cell2 == "ICM" or cell1 == "ICM"
            ):
                continue
            cell_file_1 = file_names1[methylation_type_1][cell1]
            cell_file_2 = file_names2[methylation_type_2][cell2]

            corr[i, j] = correlation(cell_file_1, cell_file_2)
            print(
                f"{cell1} ({dataset_name_1}) vs {cell2} ({dataset_name_2}): {corr[i, j]}"
            )
            ang[i, j] = correlation_to_angular_distance(corr[i, j])

    print(corr)

    # Plot correlation heatmap
    plt.figure(figsize=figsize)  # Adjust size as needed
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        xticklabels=display_names2,
        yticklabels=display_names1,
    )
    plt.title(f"Cross-dataset Correlation Matrix for {methylation_type_2}")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            RESULTS_DIR,
            f"corr_matrix_{methylation_type_2}_{dataset_name_1}_{dataset_name_2}.png",
        )
    )
    plt.close()


def correlate_everything(dataset_params=dataset_params):
    # First, create a list of all unique combinations
    combinations = []
    labels = []
    for dataset in dataset_params:
        for mtype in dataset["methylation_types"]:
            for cell_type in dataset["cell_types"]:
                if mtype == "GCA.GCC.GCT":
                    continue  # FIXME: as a hack to compute quicker
                combinations.append((dataset, mtype, cell_type))
                label = f"{dataset.name}_{mtype}_{cell_type}"
                # print(label)
                labels.append(label)

    n = len(combinations)
    # Create the full correlation matrix
    corr_matrix = np.zeros((n, n))

    # print(f"{combinations}\ncombinations above")

    # Fill the correlation matrix
    for i, (d1, m1, c1) in enumerate(combinations):
        for j, (d2, m2, c2) in enumerate(combinations):
            # Skip if ICM condition is met
            if (m1 == "GCA.GCC.GCT" or m2 == "GCA.GCC.GCT") and (
                c1 == "ICM" or c2 == "ICM"
            ):
                corr_matrix[i, j] = np.nan
                continue

            # Get the file names
            cell_file_1 = d1["filenames"][m1][c1]
            cell_file_2 = d2["filenames"][m2][c2]

            corr = correlation(cell_file_1, cell_file_2)
            corr_matrix[i, j] = corr
            print(f"Correlation = {corr}")

    return corr_matrix, labels


def calculate_angular_distance(correlation_matrix):
    """
    Calculate angular distance from correlation matrix.
    Angular distance = arccos(correlation) / π
    This scales the distance to [0,1] where:
    - 0 means perfectly correlated
    - 0.5 means uncorrelated
    - 1 means perfectly anti-correlated
    """
    # Handle any numerical instabilities in correlation values
    correlation_matrix = np.clip(correlation_matrix, -1.0, 1.0)

    # Calculate angular distance
    angular_dist = np.arccos(correlation_matrix) / np.pi

    # Handle any NaN values if present
    if np.any(np.isnan(angular_dist)):
        # Replace NaN with maximum distance (1.0)
        angular_dist = np.nan_to_num(angular_dist, nan=1.0)

    return angular_dist


def cluster_data(correlation_matrix, labels, method="average", metric="angular"):
    """
    Perform hierarchical clustering using angular distance.

    Parameters:
    - correlation_matrix: The correlation matrix
    - labels: List of labels for each row/column
    - method: Linkage method ('single', 'complete', 'average', 'weighted', 'ward')
    - metric: Distance metric ('angular' or 'correlation')

    Returns:
    - linkage_matrix: The hierarchical clustering linkage matrix
    - distance_matrix: The distance matrix used for clustering
    """

    if metric == "angular":
        distance_matrix = calculate_angular_distance(correlation_matrix)
    else:  # 'correlation'
        distance_matrix = 1 - correlation_matrix

    # Convert to condensed form (required by scipy)
    condensed_dist = squareform(distance_matrix)

    # Perform hierarchical clustering
    linkage_matrix = hierarchy.linkage(condensed_dist, method="single")

    return linkage_matrix, distance_matrix


# Example usage:
def plot_clustering(correlation_matrix, labels, method="average", metric="angular"):
    """
    Plot the hierarchical clustering dendrogram.
    """

    # Perform clustering
    linkage_matrix, distance_matrix = cluster_data(
        correlation_matrix, labels, method, metric
    )

    # Create dendrogram
    plt.figure(figsize=(15, 10))
    hierarchy.dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,  # rotates the x axis labels
        leaf_font_size=8,  # font size for the x axis labels
    )
    plt.title(
        f"Hierarchical Clustering Dendrogram ({method} linkage, {metric} distance)"
    )
    plt.xlabel("Sample")
    plt.ylabel("Distance")
    plt.tight_layout()

    return linkage_matrix, distance_matrix


# Optional: Get cluster assignments for a specific number of clusters
def get_clusters(linkage_matrix, n_clusters, labels):
    """
    Get cluster assignments for each sample.
    """
    from scipy.cluster import hierarchy

    clusters = hierarchy.fcluster(linkage_matrix, n_clusters, criterion="maxclust")
    cluster_assignments = dict(zip(labels, clusters))
    return cluster_assignments


def main():
    # Run the full analysis
    correlation_matrix, matrix_labels = correlate_everything(dataset_params)
    linkage_matrix, distance_matrix = plot_clustering(
        correlation_matrix, matrix_labels, method="average", metric="angular"
    )

    plt.figure(figsize=(15, 10))
    hierarchy.dendrogram(
        linkage_matrix,
        labels=matrix_labels,
        leaf_rotation=90,  # rotates the x axis labels
        leaf_font_size=8,  # font size for the x axis labels
    )
    plt.title(f"Hierarchical Clustering Dendrogram Methylation Data")
    plt.xlabel("Sample")
    plt.ylabel("Distance (angle)")
    plt.tight_layout()
    plt.show()
    for params in dataset_params:
        generate_sample_histograms(params)
    sample_correlations()
    sample_correlations(dataset="WGBS", methylation_type="ACG.TCG")


main()
