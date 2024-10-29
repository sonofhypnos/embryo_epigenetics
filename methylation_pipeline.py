#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import json
from datetime import datetime
import multiprocessing
import glob

cpu_count = multiprocessing.cpu_count()

cpu_count_safety_margin = max([cpu_count - 4, (cpu_count * 3) // 4 ])

#Note that it is not recommended to remove too-short sequences if the analysed FastQ file is one of a pair of paired-end files, since this confuses the sequence-by-sequence order of paired-end reads which is again required by many aligners. For paired-end files, Trim Galore! has an option --paired which runs a paired-end validation on both trimmed _1 and _2 FastQ files once the trimming has completed. This step removes entire read pairs if at least one of the two sequences became shorter than a certain threshold. If only one of the two reads is longer than the set threshold, e.g. when one read has very poor qualities throughout, this singleton read can be written out to unpaired files (see option retain_unpaired) which may be aligned in a single-end manner.

# TODO: add really simple tests for some of these steps
# TODO: Use trim-galore properly (see documentation trim-galore): "for paired-end files, Trim
# Galore! has an option --paired which runs a paired-end validation on both
# trimmed _1 and _2 FastQ files once the trimming has completed. This step
# removes entire read pairs if at least one of the two sequences became shorter
# than a certain threshold. If only one of the two reads is longer than the set
# threshold, e.g. when one read has very poor qualities throughout, this
# singleton read can be written out to unpaired files (see option
# retain_unpaired) which may be aligned in a single-end manner."
# NOTE: it seems like making a custom index for our genes is in fact not that much faster? (we do need less memory-> we can go with a higher number of cores) (seems like we were in fact still using entire chromosomes and this is why it took that long!)
# TODO: the bismark reference recommends a deduplication step for full genome
# sequencing (not for RRBS). Learn what that is useful for? (and then add it to our pipeline)
# TODO: check tradeoffs between using a more variable reference genome (leads to
# possibly non-unique reads and more computationally expensive (the non-unique
# part seems like that would be desirable, because we then know that at that
# part it was ambigous)) and using more lenient parameters for the matching with bismark
# TODO: given that my quality scores are this high, I think I want to try it
# with the alternative genome version (although perhaps first check how other
# people do this step online)
# TODO: for every shell command we use, check if we are using a parallel version
# TODO: delete fastq files after checking integrity of the trimmed files
# TODO: understand what bismark 2 bed graph does and if that is something we need
# TODO: get a deeper understanding of the types of algorithms we are using here and if we have different alternatives
# TODO: check that fq_C_to_T files that were in the root directory of the project would have needed to be somewhere else
# TODO: make sure we have enough space for intermediate files
# TODO: check if we can save space by using gziped versions of the files with bismark
# NOTE: on overall performance: it seems like the trimming is also taking quite a few minutes, but nowhere near as much time as the alignment step

class BisulfiteAnalyzer:
    def __init__(self, output_dir="methylation_analysis"):
        """
        Initialize analyzer with output directory for intermediate files
        """
        self.output_dir = output_dir
        self.checkpoint_file = os.path.join(output_dir, "checkpoint.json")
        os.makedirs(output_dir, exist_ok=True)
        
        # Load checkpoint if exists
        self.checkpoint = self.load_checkpoint()
        
    def load_checkpoint(self):
        """Load checkpoint from file if it exists"""
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {"processed_files": {}}
    
    def save_checkpoint(self, sra_file, stage):
        """Save processing progress to checkpoint file"""
        self.checkpoint["processed_files"][sra_file] = {
            "stage": stage,
            "timestamp": datetime.now().isoformat()
        }
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoint, f, indent=2)
    
    # def get_stage(self, sra_file):
    #     """Get the last completed stage for a file"""
    #     return self.checkpoint["processed_files"].get(sra_file, {}).get("stage", None)
    def get_stage(self, sra_file):
        """Return list of all completed stages for a file"""
        stages = []
        for file, info in self.checkpoint["processed_files"].items():
            if file.startswith(sra_file):  # Match files that start with the SRA ID
                stages.append(info.get("stage"))
        return stages


    def run_fastq_dump(self, sra_file):
        """Convert SRA to FASTQ format"""
        output_dir = os.path.join(self.output_dir, "fastq")
        # os.makedirs(output_dir, exist_ok=True)

        base_name = os.path.basename(sra_file)
        fastq1 = os.path.join(output_dir, f"{base_name}_1.fastq.gz")
        fastq2 = os.path.join(output_dir, f"{base_name}_2.fastq.gz")
        
        # # Skip if files exist and stage is completed
        # if os.path.exists(fastq1) and os.path.exists(fastq2) and self.get_stage(sra_file) == "fastq_dump":
        #     print(f"Skipping fastq-dump for {base_name} - already processed")
        #     return fastq1, fastq2

        # cmd = f"fastq-dump --split-files --outdir {output_dir} {sra_file}"
        # subprocess.run(cmd, shell=True, check=True)
        # self.save_checkpoint(sra_file, "fastq_dump")
        return fastq1, fastq2
    
    def run_trim_galore(self, fastq1, fastq2):
        """Trim adapters and low-quality sequences."""
        output_dir = os.path.join(self.output_dir, "trimmed")
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = os.path.basename(fastq1).replace("_1.fastq.gz", "")
        trimmed1 = os.path.join(output_dir, f"{base_name}_1_val_1.fq.gz")
        trimmed2 = os.path.join(output_dir, f"{base_name}_2_val_2.fq.gz")
        
        # Skip if files exist and stage is completed
        if os.path.exists(trimmed1) and os.path.exists(trimmed2) and "trim_galore" in self.get_stage(base_name):
            print(f"Skipping trim_galore for {base_name} - already processed")
            return trimmed1, trimmed2

        # NOTE: in the paper they note clipping 6
        # NOTE: I will see if not clipping at all is better
        #
        # cmd = f"trim_galore  --paired --length 20 --clip_R1 6 --clip_R2 6 --cores 8 --output_dir {output_dir} {fastq1} {fastq2}"
        cmd = f"trim_galore  --paired --cores 8 --output_dir {output_dir} {fastq1} {fastq2}"
        subprocess.run(cmd, shell=True, check=True)
        self.save_checkpoint(base_name, "trim_galore")
        return trimmed1, trimmed2

    def find_latest_bam(self, output_dir, base_name):
        """Find the final BAM file, excluding temporary files"""
        # Pattern for final BAM file (no .temp. in name)
        pattern = os.path.join(output_dir, f"{base_name}*_bismark_bt2_pe.bam")
        bam_files = glob.glob(pattern)

        # Filter out temporary files
        final_bams = [f for f in bam_files if '.temp.' not in f]

        if not final_bams:
            return None
        return final_bams[0]



    def run_bismark_alignment(self, trimmed1, trimmed2, reference_genome):
        output_dir = os.path.join(self.output_dir, "aligned")
        os.makedirs(output_dir, exist_ok=True)

        base_name = os.path.basename(trimmed1).replace("_1_val_1.fq", "")
        bam_file = self.find_latest_bam(output_dir, base_name)

        if bam_file and os.path.exists(bam_file) and "bismark" in self.get_stage(bam_file):
            print(f"Skipping bismark alignment for {base_name} - already processed")
            return bam_file

        # NOTE: According to the documentation (https://felixkrueger.github.io/Bismark/options/alignment/), bismark multicore roughly increases RAM requirements in a linear fashion.
        # > If system resources are plentiful this is a viable option to speed
        # > up the alignment process (we observed a near linear speed increase
        # > for up to --multicore 8 tested). However, please note that a typical
        # > Bismark run will use several cores already (Bismark itself, 2 or 4
        # > threads for Bowtie/Bowtie2, Samtools, gzip etc...) and ~10-16GB of
        # > memory per thread depending on the choice of aligner and genome.
        # > WARNING: Bismark Parallel is resource hungry! Each value of
        # > --multicore specified will effectively lead to a linear increase in
        # > compute and memory requirements, so --parallel 4 for e.g. the GRCm38
        # > mouse genome will probably use ~20 cores and eat ~48GB of RAM, but
        # > at the same time reduce the alignment time to ~25-30%. You have been
        # > warned.
        # So we will not set a higher value to not overload the memory
        # NOTE: (according to the above, 8 should already use ~90 GB, I saw it using 60GB (could have been higher when I wasn't looking))
        # FIXME: score_min might be too lenient, we should check what typical values here for bisulfite sequencing are

        cmd = f"bismark --parallel 8 --output_dir {output_dir} --genome {reference_genome} -1 {trimmed1} -2 {trimmed2} --score_min L,0,-0.6" #Using -0.6 since that is apparently a standard value
        subprocess.run(cmd, shell=True, check=True)

        # Find the generated BAM file
        bam_file = self.find_latest_bam(output_dir, base_name)
        if not bam_file:
            raise FileNotFoundError(f"BAM file not found for {base_name}")

        self.save_checkpoint(bam_file, "bismark")
        return bam_file
    
    def extract_methylation(self, bam_file, genes_of_interest, gene_coords, mapping_file):
        output_dir = os.path.join(self.output_dir, "methylation")
        os.makedirs(output_dir, exist_ok=True)

        base_name = os.path.basename(bam_file).replace("_bismark_bt2_pe.bam", "")
        methylation_file = os.path.join(output_dir, f"{base_name}_methylation.txt")

        if os.path.exists(methylation_file) and "methylation" in self.get_stage(base_name):
            print(f"Loading existing methylation data for {base_name}")
            return self.load_methylation_results(methylation_file)

        # Extract methylation data
        print(f"Extracting methylation data from {bam_file}")
        extract_cmd = f"bismark_methylation_extractor --parallel {max([cpu_count_safety_margin // 3, 1])} --comprehensive --bedGraph --output {output_dir} {bam_file}"
        subprocess.run(extract_cmd, shell=True, check=True)

        # Find the generated bedGraph file
        bedgraph_file = glob.glob(os.path.join(output_dir, f"*{base_name}*.bedGraph.gz"))[0]

        # Convert coordinates to genomic positions
        converted_bedgraph = os.path.join(output_dir, f"{base_name}_converted.bedGraph")
        converted_file = convert_methylation_coordinates(bedgraph_file, mapping_file, converted_bedgraph)

        # Use converted coordinates for intersection
        methylation_data = {}
        for gene in genes_of_interest:
            if not os.path.exists(gene_coords[gene]):
                print(f"Warning: BED file not found for {gene}: {gene_coords[gene]}")
                continue

            print(f"Analyzing methylation for {gene}")
            cmd = f"bedtools intersect -a {converted_file} -b {gene_coords[gene]}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

            if result.stderr:
                print(f"Warning for {gene}: {result.stderr}")

            methylation_data[gene] = self.parse_methylation_output(result.stdout)
            print(f"Found {methylation_data[gene]['total_sites']} sites for {gene}")

        # Save results
        with open(methylation_file, 'w') as f:
            json.dump(methylation_data, f, indent=2)

        self.save_checkpoint(base_name, "methylation")
        return methylation_data

    @staticmethod
    def parse_methylation_output(methylation_text):
        """Parse methylation data and calculate statistics"""
        total_sites = 0
        methylated_count = 0
        total_coverage = 0

        for line in methylation_text.splitlines():
            if not line.strip():
                continue

            fields = line.split('\t')
            if len(fields) >= 4:  # Ensure we have enough fields
                coverage = int(fields[3])
                methyl_percent = float(fields[4])

                total_sites += 1
                total_coverage += coverage
                methylated_count += (coverage * methyl_percent / 100)

        return {
            'total_sites': total_sites,
            'total_coverage': total_coverage,
            'methylated_sites': methylated_count,
            'methylation_percentage': (methylated_count/total_coverage*100) if total_coverage > 0 else 0
        }

    @staticmethod
    def load_methylation_results(methylation_file):
        with open(methylation_file, 'r') as f:
            return json.load(f)

def analyze_methylation(fastq_files, reference_genome, genes_of_interest, gene_coords, output_dir, input_dir="fastq"):
    """Main analysis function"""

    analyzer = BisulfiteAnalyzer(output_dir=output_dir)
    results = {}
    
    for fastq_file in fastq_files:
        base_name = os.path.basename(fastq_file)
        print(f"\nProcessing {base_name}...")


        fastq1 = fastq_file
        fastq2 = fastq1.replace("_1.fastq.gz", "_2.fastq.gz")
        assert os.path.exists(fastq1)
        print(fastq2)
        assert os.path.exists(fastq2)

        print(f"Done with fastq conversion: {fastq1}, {fastq2}")
        # Quality control and trimming (if needed)
        trimmed1, trimmed2 = analyzer.run_trim_galore(fastq1, fastq2)
        print(f"Done with Quality control and trimming")
        # FIXME: only here because alignment is currently broken:
        continue

        # Alignment (if needed)
        bam_file = analyzer.run_bismark_alignment(trimmed1, trimmed2, reference_genome)
        print(f"done with alignment")
        
        # Methylation analysis (if needed)
        results[base_name] = analyzer.extract_methylation(bam_file, genes_of_interest, gene_coords)
    
    return results

if __name__ == "__main__":
    genes_of_interest = ['NANOG', 'KLF17', 'DPPA4', 'NANOG_promoter', 'KLF17_promoter', 'DPPA4_promoter']
    
    # fastq_pattern = os.path.join(f"fastq/*_1.fastq.gz")
    # fastq_files = glob.glob(fastq_pattern)

    # TODO: add special case for sperm
    # TODO: check we didn't miss any SRR files (check if left over files under ~/sra_data/ are present as fastq files)
    fastq_files = [os.path.join(f"fastq/SRR6228411_1.fastq.gz")]

    reference_genome = '/home/ubuntu/sra_data/reference_genome/'
    gene_coords = {
        'NANOG': 'bedfiles/NANOG.bed',
        'KLF17': 'bedfiles/KLF17.bed',
        'DPPA4': 'bedfiles/DPPA4.bed',
        'NANOG_promoter': 'bedfiles/NANOG_promoter.bed',
        'KLF17_promoter': 'bedfiles/KLF17_promoter.bed',
        'DPPA4_promoter': 'bedfiles/DPPA4_promoter.bed',
    }
    version = "1.0" # Added versioning to try different parameters

    assert len(genes_of_interest) == len(gene_coords)
    assert sorted(genes_of_interest) == sorted([x for x in gene_coords.keys()])

    # Run analysis
    results = analyze_methylation(fastq_files, reference_genome, genes_of_interest, gene_coords, f"methylation_analysis_{version}")
    
    # Create summary DataFrame
    summary_df = pd.DataFrame.from_dict(results, orient='index')
    print(f"\nMethylation Analysis Summary (file:{results}):")
    print(summary_df)
