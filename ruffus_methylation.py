#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import json
from datetime import datetime
import multiprocessing
import glob
from ruffus import *
from ruffus import formatter as ruffus_formatter  # Add this line
import logging
import sys


# Pipeline parameters
PARAMS = {
    'reference_genome': '/home/ubuntu/sra_data/reference_genome/',
    'output_dir': 'methylation_analysis_2.0',
    'max_parallel_jobs': 32
}

# Create output directory if it doesn't exist
os.makedirs(PARAMS['output_dir'], exist_ok=True)

# Configure logging to both file and stdout
log_file = os.path.join(PARAMS['output_dir'], 'pipeline.log')

# Create a formatter - using logging.Formatter() instead of Formatter()
formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# Create file handler
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.INFO)
file_handler.setFormatter(formatter)

# Create console handler
console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)

# Configure root logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# TODO: add really simple tests for some of these steps
# NOTE: it seems like making a custom index for our genes is in fact not that much faster? (we do need less memory-> we can go with a higher number of cores) (seems like we were in fact still using entire chromosomes and this is why it took that long!)
# TODO: the bismark reference recommends a deduplication step for full genome
# sequencing (not for RRBS). Learn what that is useful for? (and then add it to our pipeline)
# NOTE: we do not know why, but testing on a small subset of the dataset
# revealed using non-directional mapping with bismark to be if anything better on our dataset (less bias in matching sequences and higher mapq scores)
# TODO: Do further quality control on all our sequences (so far I only did quality control on small samples)
# TODO: Add befile step
def setup_output_dirs():
    """Create output directories"""
    dirs = ['trimmed', 'aligned']
    for d in [os.path.join(PARAMS['output_dir'], subdir) for subdir in dirs]:
        os.makedirs(d, exist_ok=True)

# NOTE: this part of the pipeline hasn't been run so far. I had problems with the regex for the files
@transform("fastq/*_1.fastq.gz",
          ruffus_formatter("fastq/(?P<SID>[^/]+)_1.fastq.gz"),
          [os.path.join(PARAMS['output_dir'], "trimmed", "{SID[0]}_1_val_1.fq.gz"),
           os.path.join(PARAMS['output_dir'], "trimmed", "{SID[0]}_2_val_2.fq.gz")],
          "{SID[0]}")
def trim_reads(input_r1, output_files, sample_name):
    """Trim adapters, low-quality sequences, and polyA/polyT sequences"""
    logger.info(f"Starting trimming for sample: {sample_name}")
    
    # Setup input/output paths
    input_r2 = input_r1.replace("_1.fastq.gz", "_2.fastq.gz")
    if not os.path.exists(input_r2):
        raise ValueError(f"Missing R2 file for {input_r1}")
    
    output_dir = os.path.join(PARAMS['output_dir'], "trimmed")
    os.makedirs(output_dir, exist_ok=True)
    
    output_r1, output_r2 = output_files
    
    # Create temporary file paths
    temp_dir = os.path.join(output_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    temp1 = os.path.join(temp_dir, f"{sample_name}_temp1.fq.gz")
    temp2 = os.path.join(temp_dir, f"{sample_name}_temp2.fq.gz")
    

    # NOTE: in the paper they note clipping 6
    # NOTE: I will see if not clipping at all is better
    # NOTE: --polyA should filter the polyA sequences. This setting was reported as experimental though, so maybe look for alternatives!
    # NOTE: we might still want to save polyA sequences (like in this example:)
    # PLEASE NOTE: The poly-A trimming mode expects that sequences were both adapter and quality trimmed
                    # before looking for Poly-A tails, and it is the user's responsibility to carry out an initial round of
                    # trimming. The following sequence:

                    # 1) trim_galore file.fastq.gz
                    # 2) trim_galore --polyA file_trimmed.fq.gz
                    # 3) zcat file_trimmed_trimmed.fq.gz | grep -A 3 PolyA | grep -v ^-- > PolyA_trimmed.fastq

    # FIXME: we don't know why we have so many polyA sequences in our data! (check why)
    # cmd = f"trim_galore  --paired --length 20 --clip_R1 6 --clip_R2 6 --cores 8 --output_dir {output_dir} {fastq1} {fastq2}"
    # NOTE: we clip the sequences since we saw bias in the sequences when using fastqc (and the original paper also made this decision)
    # Initial trim_galore command
    trim_cmd = (
        f"trim_galore --paired "
        f"--clip_R1 6 --clip_R2 6 "
        f"--three_prime_clip_R1 1 --three_prime_clip_R2 1 "
        f"--output_dir {output_dir} "
        f"{input_r1} {input_r2}"
    )
    
    # Second trim_galore for polyA
    poly_cmd = (
        f"trim_galore --paired --polyA "
        f"{output_r1} {output_r2}"
    )
    
    # Cutadapt commands: (apparently polyA doesnt cut out all poly sequences)
    cut_polyT_cmd = (
        f"cutadapt -a 'T{{10}}' -A 'T{{10}}' -g 'T{{10}}' -G 'T{{10}}' -m 20 "
        f"-o {temp1} -p {temp2} {output_r1} {output_r2} "
        f"&& mv {temp1} {output_r1} "
        f"&& mv {temp2} {output_r2}"
    )
    
    cut_polyA_cmd = (
        f"cutadapt -a 'A{{10}}' -A 'A{{10}}' -g 'A{{10}}' -G 'A{{10}}' -m 20 "
        f"-o {temp1} -p {temp2} {output_r1} {output_r2} "
        f"&& mv {temp1} {output_r1} "
        f"&& mv {temp2} {output_r2}"
    )
    
    try:
        # Execute commands in sequence
        logger.info("Running initial trim_galore")
        subprocess.run(trim_cmd, shell=True, check=True)
        
        logger.info("Running polyA trim_galore")
        subprocess.run(poly_cmd, shell=True, check=True)
        
        logger.info("Running polyT cutadapt")
        subprocess.run(cut_polyT_cmd, shell=True, check=True)
        
        logger.info("Running polyA cutadapt")
        subprocess.run(cut_polyA_cmd, shell=True, check=True)
        
        # Clean up temp directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
            
        logger.info(f"Completed trimming for sample: {sample_name}")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error in trim_reads for {sample_name}: {str(e)}")
        # Clean up any partial outputs
        for f in [output_r1, output_r2, temp1, temp2]:
            if os.path.exists(f):
                os.remove(f)
        raise
    except Exception as e:
        logger.error(f"Unexpected error in trim_reads for {sample_name}: {str(e)}")
        raise

@follows(trim_reads)
@transform(trim_reads,
          ruffus_formatter(".+/(?P<basename>.+)_1_val_1\.fq\.gz$"),
          os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_bismark_bt2_pe.bam"))
# @transform(trim_reads,
#           ruffus_formatter("_1_val_1.fq.gz$"),
#           os.path.join(PARAMS['output_dir'], "aligned", "{basename[0]}_bismark_bt2_pe.bam"))
def align_bismark(input_files, output_file):
    """Align reads using Bismark"""
    trimmed1, trimmed2 = input_files
    output_dir = os.path.join(PARAMS['output_dir'], "aligned")

    # We chose the score for score_min, because this value still gave sequences where the majority 75% of errors was due to C->T or A->G errors
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
            multiprocess=PARAMS["max_parallel_jobs"],
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
