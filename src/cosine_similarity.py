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

PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
WGBS_DATA_DIR = f"{PROJECT_DIR}data/"
RRBS_DATA_DIR = f"{PROJECT_DIR}rrbs_data/"
CACHE_DIR = f"{PROJECT_DIR}tmp/"
RESULTS_DIR = f"{PROJECT_DIR}results/"

cache = Cache(CACHE_DIR)

# TODO: Add other methylation types from RRBS data


CELL_TYPES_RRBS = [
    "MII_Oocyte",
    "MII_Oocyte+scRRBS",
    "Sperm",
    "Sperm+scRRBS",
    "Zygote",
    "2-cell",
    "4-cell",
    "8-cell",
    "Morula",
    "ICM",
    "ICM+WGBS",
    "TE",
    "Liver+WGBS",
    "2nd_PB",
    "1st_PB",
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
]

methylation_types_wgbs = ["ACG.TCG", "GCA.GCC.GCT"]
methylation_types_rrbs = ["CpG"]

# TODO: handle scRRBS data once I know what that is
# TODO: handle files where we have RRBS and WGBS like ICM?
cell_type_to_seq_type = {
    "MII_Oocyte": "RRBS",
    "MII_Oocyte+scRRBS": "scRRBS",
    "Sperm": "RRBS",
    "Sperm+scRRBS": "scRRBS",
    "Zygote": "RRBS",
    "2-cell": "RRBS",
    "4-cell": "RRBS",
    "8-cell": "RRBS",
    "Morula": "RRBS",
    "ICM": "RRBS",  # WGBS also possible
    "ICM+WGBS": "WGBS",  # WGBS also possible
    "TE": "RRBS",
    "Liver+WGBS": "WGBS",
    "2nd_PB": "RRBS",
    "1st_PB": "RRBS",
    "PN+scRRBS": "scRRBS",
}


wgbs_filenames = {
    methylation_type: {
        cell_type: f"{WGBS_DATA_DIR}{cell_type}_{methylation_type}.wig"
        for cell_type in CELL_TYPES_WGBS
    }
    for methylation_type in methylation_types_wgbs
}

rrbs_filenames = {
    methylation_type: {
        cell_type: f"{RRBS_DATA_DIR}{cell_type.split('+')[0]}_{methylation_type}_{cell_type_to_seq_type[cell_type]}.wig"
        for cell_type in CELL_TYPES_RRBS
    }
    for methylation_type in methylation_types_rrbs
}


# NOTE: I don't have the storage for this one (nor the data bandwidth)! (50GB)
# wgbs_sample_filenames = {
#     methylation_type: {
#         cell_type: f"{WGBS_DATA_DIR}{cell_type}_{methylation_type}.bw"
#         for cell_type in CELL_TYPES_WGBS
#     }
#     for methylation_type in methylation_types_wgbs
# }

rrbs_sample_filenames = {
    methylation_type: {
        cell_type: glob.glob(
            f"{RRBS_DATA_DIR}GS*{cell_type_to_seq_type[cell_type]}_{cell_type.split('+')[0]}*_methylation_calling_{methylation_type}.bw"
        )
        for cell_type in CELL_TYPES_RRBS
    }
    for methylation_type in methylation_types_rrbs
}

# print(wgbs_filenames)
# print(rrbs_filenames)
# print(rrbs_sample_filenames)

dataset_params = {
    "rrbs": {
        "methylation_types": methylation_types_rrbs,
        "filenames": rrbs_filenames,
        "cell_types": CELL_TYPES_RRBS,
        "sample_filenames": rrbs_sample_filenames,
    },
    "wgbs": {
        "methylation_types": methylation_types_wgbs,
        "filenames": wgbs_filenames,
        "cell_types": CELL_TYPES_WGBS,
    },
}


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


@cache.memoize()
def correlation(f1, f2, capture_output=True):
    try:
        # results = subprocess.run(
        #     ["zsh", "-c", f"{PROJECT_DIR}src/wig_correlation {f1} {f2}"],
        #     capture_output=capture_output,
        #     text=True,
        #     check=False,
        # )

        results = run_conda_command(f"wigCorrelate {f1} {f2}")
        corr = float(results.stdout.split()[-1])
    except TypeError as e:
        print(results)
        print(e)
    except Exception as e:
        print(e)
        return 0

    # print(result)

    return corr


def sample_correlations(dataset="RRBS"):
    params = dataset_params[dataset]
    methylation_type = "CpG"
    for cell_type in params["cell_types"]:
        correlations = run_conda_command(
            f"wigCorrelate {' '.join(params['sample_filenames'][methylation_type][cell_type])}"
        )

        samples = [sample.split("\t") for sample in correlations.split("\n")][
            :-1
        ]  # removing empty end

        files = list(set([file for sample in samples for file in sample[:2]]))
        indexes = [name.split("_")[0] for name in files]
        length = len(files)
        corr = np.zeros((length, length))
        for sample in samples:
            i = files.index(sample[0])
            j = files.index(sample[1])
            val = float(sample[2])
            corr[i, j] = val
            corr[j, i] = val

        figsize = get_figure_size(length, length)

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
            f"{RESULTS_DIR}corr_matrix_{methylation_type}_{dataset}_{cell_type}.png"
        )


sample_correlations()

for dataset in dataset_params:
    for methylation_type in dataset_params[dataset]["methylation_types"]:
        cell_types = dataset_params[dataset]["cell_types"]
        file_names = dataset_params[dataset]["filenames"]

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

                print(cell_file_1)
                print(cell_file_2)
                # pearson = run_conda_command(
                #     f"wiggletools pearson trim {cell_file_1} {cell_file_2}"
                # )[1].split("\n")[0]

                # corr[i, j] = pearson
                corr[i, j] = correlation(cell_file_1, cell_file_2)
                print(cell1, cell2, corr[i, j])

                ang[i, j] = correlation_to_angular_distance(corr[i, j])
        print(corr)

        plt.figure(figsize=figsize)
        sns.heatmap(
            corr,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            xticklabels=cell_types,
            yticklabels=cell_types,
        )
        plt.title(f"Correlation Matrix for {methylation_type} for {dataset} dataset")
        plt.savefig(f"{RESULTS_DIR}correlation_matrix_{methylation_type}_{dataset}.png")
        # plt.show()

        plt.figure(figsize=figsize)
        sns.heatmap(
            ang,
            annot=True,
            fmt=".2f",
            cmap="coolwarm",
            xticklabels=cell_types,
            yticklabels=cell_types,
        )
        plt.title(f"Angle Matrix for {methylation_type} for {dataset} dataset")
        plt.savefig(f"{RESULTS_DIR}angle_matrix_{methylation_type}_{dataset}.png")
        # plt.show()


# For cell types that exist in both datasets, create display names
def get_display_name(cell_type, dataset):
    return f"{cell_type}_{dataset}"


dataset_name_1 = "wgbs"
dataset_name_2 = "rrbs"
# Loop through methylation types that exist in both datasets
# common_methylation_types = set(datasets[dataset_name_1]["methylation_types"]) & set(
#     datasets[dataset_name_2]["methylation_types"]
# )

cell_types1 = dataset_params[dataset_name_1]["cell_types"]
cell_types2 = dataset_params[dataset_name_2]["cell_types"]
file_names1 = dataset_params[dataset_name_1]["filenames"]
file_names2 = dataset_params[dataset_name_2]["filenames"]
methylation_type_1 = "ACG.TCG"
methylation_type_2 = "CpG"

# Create matrix of size (len(cell_types1) x len(cell_types2))
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
        if methylation_type_1 == "GCA.GCC.GCT" and (cell2 == "ICM" or cell1 == "ICM"):
            continue
        cell_file_1 = file_names1[methylation_type_1][cell1]
        cell_file_2 = file_names2[methylation_type_2][cell2]
        print(cell_file_1)
        print(cell_file_2)

        corr[i, j] = correlation(cell_file_1, cell_file_2)
        print(f"{cell1} ({dataset_name_1}) vs {cell2} ({dataset_name_2}): {corr[i, j]}")
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
    f"{RESULTS_DIR}corr_matrix_{methylation_type_2}_{dataset_name_1}_{dataset_name_2}.png"
)

# Plot angle heatmap
plt.figure(figsize=figsize)
sns.heatmap(
    ang,
    annot=True,
    fmt=".2f",
    cmap="coolwarm",
    xticklabels=display_names2,
    yticklabels=display_names1,
)
plt.title(f"Cross-dataset Angle Matrix for {methylation_type_2}")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(
    f"{RESULTS_DIR}ang_matrix_{methylation_type_2}_{dataset_name_1}_{dataset_name_2}.png"
)


# def wiggle_product(v1, v2):
#     run_conda_command("wiggletools product {v1} {v2}")
