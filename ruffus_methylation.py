#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import json
from datetime import datetime
import multiprocessing
import glob
from ruffus import *
from ruffus import formatter as ruffus_formatter
import sys
import argparse
import logging
import shutil

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='RNA-seq analysis pipeline')
    parser.add_argument('--test', action='store_true',
                       help='Run pipeline on test files only (files starting with "test" in /fastq/)')
    parser.add_argument('--test2', action='store_true',
                       help='Run pipeline on test files only (files starting with "test" in /fastq/)')
    parser.add_argument('--log', default='pipeline.log',
                       help='Log file path (default: pipeline.log)')
    parser.add_argument('--verbose', type=int, default=3,
                       help='Verbosity level (default: 3)')
    return parser.parse_args()

OPTIONS = parse_arguments()

(logger,mutex) = cmdline.setup_logging(
    __name__,
    OPTIONS.log,
    OPTIONS.verbose,
)

# Pipeline parameters
PARAMS = {
    'reference_genome': '/home/ubuntu/sra_data/reference_genome/',
    'output_dir': 'methylation_analysis_2.0',
    'max_parallel_jobs': 8
}

def setup_output_dirs():
    """Create output directories"""
    dirs = ['trimmed', 'trimmed/initial', 'trimmed/poly', 'trimmed/polyT', 'trimmed/polyA', 'aligned']
    for d in [os.path.join(PARAMS['output_dir'], subdir) for subdir in dirs]:
        os.makedirs(d, exist_ok=True)

def get_input_files():
    """Get input files based on test flag"""
    if OPTIONS.test:
        return [(glob.glob("fastq/test*_1.fastq.gz")[0],glob.glob("fastq/test*_2.fastq.gz")[0])]

    if OPTIONS.test2:
        return [("fastq/SRR5720828_1.fastq.gz", "fastq/SRR5720828_2.fastq.gz")]
    r1_files = glob.glob("fastq/*_1.fastq.gz")
    # For each R1 file, identify the corresponding R2 file by replacing "_1.fastq.gz" with "_2.fastq.gz"
    pairs = [(r1, r1.replace("_1.fastq.gz", "_2.fastq.gz")) for r1 in r1_files if os.path.exists(r1.replace("_1.fastq.gz", "_2.fastq.gz"))]
    return pairs


input_files = get_input_files()

@transform(input_files,
          ruffus_formatter("fastq/(?P<SID>[^/]+)_[12].fastq.gz"),
          [os.path.join(PARAMS['output_dir'], "trimmed/clean", "{SID[0]}_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed/clean", "{SID[0]}_2.fq.gz")],
          "{SID[0]}")
def initial_trim(input_files, output_files, sample_name):
    """Combined trimming and cleanup step"""
    logger.info(f"Starting trimming and cleanup for sample: {sample_name}")

    # Setup input files
    input_r1 = input_files[0]
    input_r2 = input_files[1]
    if not os.path.exists(input_r2):
        raise ValueError(f"Missing R2 file for {input_r1}")

    # Create temporary directory for intermediate files
    temp_dir = os.path.join(PARAMS['output_dir'], "trimmed/initial")
    final_dir = os.path.dirname(output_files[0])
    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(final_dir, exist_ok=True)

    # Calculate intermediate filenames (what trim_galore will create)
    temp_output_r1 = os.path.join(temp_dir, f"{sample_name}_1_val_1.fq.gz")
    temp_output_r2 = os.path.join(temp_dir, f"{sample_name}_2_val_2.fq.gz")

    # Run trim_galore
    trim_cmd = (
        f"trim_galore --paired "
        f"--clip_R1 6 --clip_R2 6 "
        f"--three_prime_clip_R1 1 --three_prime_clip_R2 1 "
        f"--output_dir {temp_dir} "
        f"{input_r1} {input_r2}"
    )

    try:
        # Run trimming
        subprocess.run(trim_cmd, shell=True, check=True)
        logger.info(f"Completed trimming for sample: {sample_name}")

        # Move files to final location with clean names
        shutil.move(temp_output_r1, output_files[0])
        shutil.move(temp_output_r2, output_files[1])
        logger.info(f"Completed filename cleanup for sample: {sample_name}")

        # Clean up temporary directory if it's empty
        if not os.listdir(temp_dir):
            os.rmdir(temp_dir)

    except subprocess.CalledProcessError as e:
        logger.error(f"Error in trimming for {sample_name}: {str(e)}")
        raise
    except Exception as e:
        logger.error(f"Error in cleanup for {sample_name}: {str(e)}")
        raise

@jobs_limit(8)
@follows(initial_trim)
@transform(initial_trim,
          formatter(".+/([^/_]+)_.*\.fq\.gz"),  # Captures everything before first underscore
          [
           # os.path.join(PARAMS['output_dir'], "aligned", "{1[0]}_bismark_bt2_pe.bam"),
           # os.path.join(PARAMS['output_dir'], "aligned", "{1[0]}_2_bismark_bt2_pe.bam"),
           os.path.join(PARAMS['output_dir'], "aligned", "{1[0]}_1.fq.gz_unmapped_reads_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "aligned", "{1[0]}_2.fq.gz_unmapped_reads_2.fq.gz")])
def align_bismark(input_files, output_files):
    """Align reads using Bismark"""
     # TODO: add really simple tests for some of these steps
     # NOTE: it seems like making a custom index for our genes is in fact not that much faster?
     # (we do need less memory-> we can go with a higher number of cores)
     # (seems like we were in fact still using entire chromosomes and this is why it took that long!)
     # TODO: the bismark reference recommends a deduplication step for full genome
     # sequencing (not for RRBS). Learn what that is useful for? (and then add it to our pipeline)
     # NOTE: we do not know why, but testing on a small subset of the dataset
     # revealed using non-directional mapping with bismark to be if anything better on our dataset
     # (less bias in matching sequences and higher mapq scores)
     # TODO: Do further quality control on all our sequences (so far I only did quality control on small samples)
     # TODO: Add bedfile step
     # TODO: perform deduplicate_bismark (to remove duplicates from pcr
     # amplification for whole genome sequencing) I think there weren't too many
     # duplicates, so this shouldn't be bad to skip.

    trimmed1, trimmed2 = input_files
    #output_bam,
    output_unmapped1, output_unmapped2 = output_files
    output_dir = os.path.join(PARAMS['output_dir'], "aligned")
    sample_name = os.path.basename(trimmed1).split('_')[0]

    logger.info(f"Starting alignment for sample: {sample_name}")

    cmd = (
        f"bismark "
        f"--output_dir {output_dir} "
        f"--genome {PARAMS['reference_genome']} "
        f"-1 {trimmed1} -2 {trimmed2} "
        f"--non_directional --unmapped"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        # Move the output BAM file to the expected location
        # bam_file = os.path.join(output_dir,
        #                        os.path.basename(trimmed1).replace('_1.fq.gz', '_bismark_bt2_pe.bam'))
        # if os.path.exists(bam_file):
        #     os.rename(bam_file, output_bam)

        # # Ensure unmapped files are in the expected location
        # for i, unmapped_file in enumerate([output_unmapped1, output_unmapped2], 1):
        #     src = os.path.join(output_dir, f"{os.path.basename(trimmed1).replace('_1.fq.gz', '')}_{i}.fastq_unmapped_reads_{i}.fq.gz")
        #     if os.path.exists(src):
        #         os.rename(src, unmapped_file)

        # logger.info(f"Completed alignment for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in align_bismark for {sample_name}: {str(e)}")
        raise

@jobs_limit(8)
@follows(align_bismark)
@transform(align_bismark,
           formatter(),
          [os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_unmapped1_bt2.bam"),
           os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_unmapped2_bt2.bam")])
def align_unmapped(input_files, output_files):
    """Align unmapped reads individually using Bismark"""
    unmapped1, unmapped2 = input_files
    output_bam1, output_bam2 = output_files
    sample_name = os.path.basename(unmapped1).split('_')[0]
    output_dir = os.path.join(PARAMS['output_dir'], "aligned")

    logger.info(f"Starting alignment of unmapped reads for sample: {sample_name}")

    # Process each unmapped file separately
    for unmapped_file, output_bam in [(unmapped1, output_bam1), (unmapped2, output_bam2)]:
        cmd = (
            f"bismark "
            f"--output_dir {output_dir} "
            f"--genome {PARAMS['reference_genome']} "
            f"{unmapped_file} "
            f"--non_directional"
        )

        subprocess.run(cmd, shell=True, check=True)

        # try:
        #     # Move the output BAM file to the expected location
        #     bam_file = os.path.join(output_dir,
        #                           os.path.basename(unmapped_file).replace('_unmapped_reads_1.fq.gz', '_bismark_bt2.bam')
        #                           if '_1.' in unmapped_file else
        #                           os.path.basename(unmapped_file).replace('_unmapped_reads_2.fq.gz', '_bismark_bt2.bam'))
        #     if os.path.exists(bam_file):
        #         os.rename(bam_file, output_bam)

        # except subprocess.CalledProcessError as e:
        #     logger.error(f"Error in align_unmapped for {sample_name}: {str(e)}")
        #     raise

    logger.info(f"Completed alignment of unmapped reads for sample: {sample_name}")

def main():
    """Main pipeline execution"""
    setup_output_dirs()

    pipeline_options = {
        'logger': logger,
    }

    try:
        pipeline_run(
            multiprocess=PARAMS['max_parallel_jobs'],
            **pipeline_options
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        print("failed!")
        raise

if __name__ == "__main__":
    main()
