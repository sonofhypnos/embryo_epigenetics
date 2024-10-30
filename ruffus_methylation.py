#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import json
from datetime import datetime
import multiprocessing
import glob
from ruffus import *
from ruffus.combinatorics import product
import logging

# Configure logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Pipeline parameters
PARAMS = {
    'reference_genome': '/home/ubuntu/sra_data/reference_genome/',
    'output_dir': 'methylation_analysis',
    'cores_per_job': 8,
    'max_parallel_jobs': 16
}

def setup_output_dirs():
    """Create output directories"""
    dirs = ['trimmed', 'aligned']
    for d in [os.path.join(PARAMS['output_dir'], subdir) for subdir in dirs]:
        os.makedirs(d, exist_ok=True)

@transform("fastq/*_1.fastq.gz",
           regex(r"fastq\/([^/]+)_1\.fastq\.gz"),
           [os.path.join(PARAMS['output_dir'], "trimmed", r"\1_1_val_1.fq.gz"),
            os.path.join(PARAMS['output_dir'], "trimmed", r"\1_2_val_2.fq.gz")],
           r"\1")
def trim_reads(input_r1, output_files, sample_name):
    """Trim adapters and low-quality sequences"""
    print(f"INPUT1: {input_r1}")
    return
    input_r2 = input_r1.replace("_1.fastq.gz", "_2.fastq.gz")
    output_dir = os.path.join(PARAMS['output_dir'], "trimmed")

    if not os.path.exists(input_r2):
        raise ValueError(f"Missing R2 file for {input_r1}")

    # Define output file paths
    output_r1, output_r2 = output_files

    # Run trim_galore with polyA trimming and output to output_dir
    cmd = (
        f"trim_galore --paired --cores {PARAMS['cores_per_job']} --clip_R1 6 --clip_R2 6 "
        f"--three_prime_clip_R1 1 --three_prime_clip_R2 1 "
        f"--output_dir {output_dir} --polyA "
        f"{input_r1} {input_r2}"
    )

    # Execute the command
    try:
        subprocess.run(cmd, shell=True, check=True)
        # Move the output files to match the expected output filenames
        os.rename(os.path.join(output_dir, f"{sample_name}_1_val_1.fq.gz"), output_r1)
        os.rename(os.path.join(output_dir, f"{sample_name}_2_val_2.fq.gz"), output_r2)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in trim_reads for {sample_name}: {str(e)}")
        raise

# @transform(trim_reads,
#            regex(r".+/(.+)_1_val_1\.fq\.gz"),
#            add_inputs(r".+/\1_2_val_2.fq.gz"),
#            os.path.join(PARAMS['output_dir'], "aligned", r"\1_bismark_bt2_pe.bam"))
@follows(trim_reads)
@transform(trim_reads,
          formatter("_1_val_1.fq.gz$"),
          os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_bismark_bt2_pe.bam"))
def align_bismark(input_files, output_file):
    """Align reads using Bismark"""
    trimmed1, trimmed2 = input_files
    output_dir = os.path.join(PARAMS['output_dir'], "aligned")

    cmd = (
        f"bismark --parallel {PARAMS['cores_per_job']} "
        f"--output_dir {output_dir} "
        f"--genome {PARAMS['reference_genome']} "
        f"-1 {trimmed1} -2 {trimmed2} "
        f"--score_min L,0,-0.6 --non_directional"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        # Move the output BAM file to the expected location
        bam_file = os.path.join(output_dir, os.path.basename(trimmed1).replace('_1_val_1.fq.gz', '_bismark_bt2_pe.bam'))
        if os.path.exists(bam_file):
            os.rename(bam_file, output_file)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in align_bismark for {trimmed1}: {str(e)}")
        raise

def main():
    """Main pipeline execution"""
    setup_output_dirs()

    # Pipeline options
    pipeline_options = {
        # 'checksum_file': os.path.join(PARAMS['output_dir'], 'ruffus_checksum.txt'),
        # 'log_file': os.path.join(PARAMS['output_dir'], 'pipeline.log'),
        'verbose': 3,
    }

    try:
        # Run the pipeline
        pipeline_run(
            target_tasks=[align_bismark],
            **pipeline_options
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        raise

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        raise
