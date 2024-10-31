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
# import logging
import sys
import shutil
import argparse

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='RNA-seq analysis pipeline')
    parser.add_argument('--test', action='store_true',
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
    'max_parallel_jobs': 32
}

def setup_output_dirs():
    """Create output directories"""
    dirs = ['trimmed', 'trimmed/initial', 'trimmed/poly', 'trimmed/polyT', 'trimmed/polyA', 'aligned']
    for d in [os.path.join(PARAMS['output_dir'], subdir) for subdir in dirs]:
        os.makedirs(d, exist_ok=True)

def get_input_files():
    """Get input files based on test flag"""
    if OPTIONS.test:
        return glob.glob("fastq/test*_1.fastq.gz")
    return glob.glob("fastq/*_1.fastq.gz")


# NOTE: this part of the pipeline hasn't been run so far. I had problems with the regex for the files
# @transform(glob.glob("fastq/*_1.fastq.gz") if not (OPTIONS and OPTIONS.test)
#           else glob.glob("fastq/test*_1.fastq.gz"),
#           ruffus_formatter("fastq/(?P<SID>[^/]+)_1.fastq.gz"),
#           [os.path.join(PARAMS['output_dir'], "trimmed/initial", "{SID[0]}_1_val_1.fq.gz"),
#            os.path.join(PARAMS['output_dir'], "trimmed/initial", "{SID[0]}_2_val_2.fq.gz")],
#            "{SID[0]}")
@transform(get_input_files,
          regex(r"fastq/(.+)_1.fastq.gz"),
          [os.path.join(PARAMS['output_dir'], "trimmed/initial", r"\1_1_val_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed/initial", r"\1_2_val_2.fq.gz")],
          r"\1")
def initial_trim(input_r1, output_files, sample_name):
    """Initial trimming step"""
    logger.info(f"Starting initial trimming for sample: {sample_name}")

    input_r2 = input_r1.replace("_1.fastq.gz", "_2.fastq.gz")
    if not os.path.exists(input_r2):
        raise ValueError(f"Missing R2 file for {input_r1}")

    output_dir = os.path.join(PARAMS['output_dir'], "trimmed/initial")

    # NOTE: we clip the sequences since we saw bias in the sequences when using fastqc
    # (and the original paper also made this decision)
    trim_cmd = (
        f"trim_galore --paired "
        f"--clip_R1 6 --clip_R2 6 "
        f"--three_prime_clip_R1 1 --three_prime_clip_R2 1 "
        f"--output_dir {output_dir} "
        f"{input_r1} {input_r2}"
    )

    try:
        subprocess.run(trim_cmd, shell=True, check=True)
        logger.info(f"Completed initial trimming for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in initial_trim for {sample_name}: {str(e)}")
        raise

@follows(initial_trim)
@transform(initial_trim,
          formatter(),
          [os.path.join(PARAMS['output_dir'], "trimmed/poly", "{basename[0]}_1_val_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed/poly", "{basename[0]}_2_val_2.fq.gz")])
def poly_trim(input_files, output_files):
    """PolyA trimming step"""
    # NOTE: we might still want to save polyA sequences (like in this example:)
    # PLEASE NOTE: The poly-A trimming mode expects that sequences were both adapter and quality trimmed
    # before looking for Poly-A tails, and it is the user's responsibility to carry out an initial round of
    # trimming.

    # FIXME: we don't know why we have so many polyA sequences in our data! (check why)

    input_r1, input_r2 = input_files
    output_r1, output_r2 = output_files
    sample_name = os.path.basename(input_r1).split('_')[0]

    logger.info(f"Starting polyA trimming for sample: {sample_name}")

    poly_cmd = (
        f"trim_galore --paired --polyA "
        f"--output_dir {os.path.dirname(output_r1)} "
        f"{input_r1} {input_r2}"
    )

    try:
        subprocess.run(poly_cmd, shell=True, check=True)
        logger.info(f"Completed polyA trimming for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in poly_trim for {sample_name}: {str(e)}")
        raise

@follows(poly_trim)
@transform(poly_trim,
          formatter(),
          [os.path.join(PARAMS['output_dir'], "trimmed/polyT", "{basename[0]}_1_val_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed/polyT", "{basename[0]}_2_val_2.fq.gz")])
def polyT_trim(input_files, output_files):
    """PolyT trimming step using cutadapt"""
    input_r1, input_r2 = input_files
    output_r1, output_r2 = output_files
    sample_name = os.path.basename(input_r1).split('_')[0]

    logger.info(f"Starting polyT trimming for sample: {sample_name}")

    cmd = (
        f"cutadapt -a 'T{{10}}' -A 'T{{10}}' -g 'T{{10}}' -G 'T{{10}}' -m 20 "
        f"-o {output_r1} -p {output_r2} {input_r1} {input_r2}"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f"Completed polyT trimming for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in polyT_trim for {sample_name}: {str(e)}")
        raise

@follows(polyT_trim)
@transform(polyT_trim,
          formatter(),
          [os.path.join(PARAMS['output_dir'], "trimmed/polyA", "{basename[0]}_1_val_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed/polyA", "{basename[0]}_2_val_2.fq.gz")])
def polyA_trim(input_files, output_files):
    """PolyA trimming step using cutadapt"""
    input_r1, input_r2 = input_files
    output_r1, output_r2 = output_files
    sample_name = os.path.basename(input_r1).split('_')[0]

    logger.info(f"Starting polyA trimming for sample: {sample_name}")

    cmd = (
        f"cutadapt -a 'A{{10}}' -A 'A{{10}}' -g 'A{{10}}' -G 'A{{10}}' -m 20 "
        f"-o {output_r1} -p {output_r2} {input_r1} {input_r2}"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        logger.info(f"Completed polyA trimming for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in polyA_trim for {sample_name}: {str(e)}")
        raise

@follows(polyA_trim)
@transform(polyA_trim,
          formatter(),
          os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_bismark_bt2_pe.bam"))
def align_bismark(input_files, output_file):
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
    # TODO: Add befile step

    trimmed1, trimmed2 = input_files
    output_dir = os.path.join(PARAMS['output_dir'], "aligned")
    sample_name = os.path.basename(trimmed1).split('_')[0]

    logger.info(f"Starting alignment for sample: {sample_name}")

    # We chose the score for score_min, because this value still gave sequences where the majority 75%
    # of errors was due to C->T or A->G errors
    cmd = (
        f"bismark "
        f"--output_dir {output_dir} "
        f"--genome {PARAMS['reference_genome']} "
        f"-1 {trimmed1} -2 {trimmed2} "
        f"--score_min L,-0.6,-1.0 --non_directional"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        # Move the output BAM file to the expected location
        bam_file = os.path.join(output_dir,
                               os.path.basename(trimmed1).replace('_1_val_1.fq.gz', '_bismark_bt2_pe.bam'))
        if os.path.exists(bam_file):
            os.rename(bam_file, output_file)

        logger.info(f"Completed alignment for sample: {sample_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in align_bismark for {sample_name}: {str(e)}")
        raise

def main():
    """Main pipeline execution"""
    setup_output_dirs()

    pipeline_options = {
        'logger': logger,
    }

    try:
        # Run the pipeline
        pipeline_run(
            target_tasks=[align_bismark],
            # multiprocess=PARAMS["max_parallel_jobs"], #NOTE: not sure this is needed
            # **pipeline_options
        )
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        print("failed!")
        raise

if __name__ == "__main__":
    # try:
    main()
    # except Exception as e:
    #     logging.error(f"Pipeline execution failed: {str(e)}")
    #     raise
