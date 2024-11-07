import subprocess
import os
from itertools import product
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


PROJECT_DIR = "/home/tassilo/repos/embryo_epigenetics/"
DATA_DIR = f"{PROJECT_DIR}data/"

CELL_TYPES_RRBS = [
    "Sperm",
    "Zygote",
    "Oocyte",
    "2-cell",
    "4-cell",
    "8-cell",
    "Morula",
    "ICM",
    "TE",
]

CELL_TYPES_WGBS = [
    "Sperm",
    "Zygote",
    "Oocyte",
    "2cell",
    "4cell",
    "8cell",
    "Morula",
    "ICM",
    "TE",
]

methylation_types_wgbs = ["ACG.TCG", "GCA.GCC.GCT"]

wgbs_filenames = {
    methylation_type: {
        cell_type: f"{DATA_DIR}{cell_type}_{methylation_type}.bw"
        for cell_type in CELL_TYPES_WGBS
    }
    for methylation_type in methylation_types_wgbs
}
print(wgbs_filenames)


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

        if capture_output:
            return result.returncode, result.stdout, result.stderr
        return result

    except subprocess.SubprocessError as e:
        if capture_output:
            return 1, "", str(e)
        raise


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
    if r < -1 or r > 1:
        raise ValueError("Correlation coefficient must be between -1 and 1")

    return np.arccos(r) / np.pi


def cosine_similarity(f1, f2):
    sum_1_sqrt = float(run_conda_command(f"wiggletools meanI mult {f1} {f1}")[1]) ** 0.5
    sum_2_sqrt = float(run_conda_command(f"wiggletools meanI mult {f2} {f2}")[1]) ** 0.5
    sum_3 = float(run_conda_command(f"wiggletools meanI mult {f1} {f2}")[1])
    return sum_3 / (sum_1_sqrt * sum_2_sqrt)


for methylation_type in methylation_types_wgbs:
    length = len(CELL_TYPES_WGBS)
    corr = np.zeros((length, length))
    ang = np.zeros((length, length))

    for i, cell1 in enumerate(CELL_TYPES_WGBS):
        for j, cell2 in enumerate(CELL_TYPES_WGBS):
            cell_file_1 = wgbs_filenames[methylation_type][cell1]
            cell_file_2 = wgbs_filenames[methylation_type][cell2]
            # pearson = run_conda_command(
            #     f"wiggletools pearson trim {cell_file_1} {cell_file_2}"
            # )[1].split("\n")[0]

            # corr[i, j] = pearson
            corr[i, j] = cosine_similarity(cell_file_1, cell_file_2)

            ang[i, j] = correlation_to_angular_distance(corr[i, j])
    print(corr)
    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        xticklabels=CELL_TYPES_WGBS,
        yticklabels=CELL_TYPES_WGBS,
    )
    plt.title("Correlation Matrix for {methylation_type}")
    plt.show()

    sns.heatmap(
        ang,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        xticklabels=CELL_TYPES_WGBS,
        yticklabels=CELL_TYPES_WGBS,
    )
    plt.title("Angle Matrix for {methylation_type}")
    plt.show()


# def wiggle_product(v1, v2):
#     run_conda_command("wiggletools product {v1} {v2}")
