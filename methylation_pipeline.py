#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
import json
from datetime import datetime
import multiprocessing

cpu_count = multiprocessing.cpu_count()

cpu_count_safety_margin = max([cpu_count - 4, (cpu_count * 3) // 4 ])

# TODO: for every shell command we use, check if we are using a parallel version
# TODO: delete fastq files after checking integrity of the trimmed files
# TODO: understand what bismark 2 bed graph does and if that is something we need
# TODO: get a deeper understanding of the types of algorithms we are using here and if we have different alternatives

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
    
    def get_stage(self, sra_file):
        """Get the last completed stage for a file"""
        return self.checkpoint["processed_files"].get(sra_file, {}).get("stage", None)

    def run_fastq_dump(self, sra_file):
        """Convert SRA to FASTQ format"""
        output_dir = os.path.join(self.output_dir, "fastq")
        # os.makedirs(output_dir, exist_ok=True)

        base_name = os.path.basename(sra_file)
        fastq1 = os.path.join(output_dir, f"{base_name}_1.fastq")
        fastq2 = os.path.join(output_dir, f"{base_name}_2.fastq")
        
        # # Skip if files exist and stage is completed
        # if os.path.exists(fastq1) and os.path.exists(fastq2) and self.get_stage(sra_file) == "fastq_dump":
        #     print(f"Skipping fastq-dump for {base_name} - already processed")
        #     return fastq1, fastq2

        # cmd = f"fastq-dump --split-files --outdir {output_dir} {sra_file}"
        # subprocess.run(cmd, shell=True, check=True)
        # self.save_checkpoint(sra_file, "fastq_dump")
        return fastq1, fastq2
    
    def run_trim_galore(self, fastq1, fastq2):
        """Trim adapters and low-quality sequences"""
        output_dir = os.path.join(self.output_dir, "trimmed")
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = os.path.basename(fastq1).replace("_1.fastq", "")
        trimmed1 = os.path.join(output_dir, f"{base_name}_1_val_1.fq")
        trimmed2 = os.path.join(output_dir, f"{base_name}_2_val_2.fq")
        
        # Skip if files exist and stage is completed
        if os.path.exists(trimmed1) and os.path.exists(trimmed2) and self.get_stage(base_name) == "trim_galore":
            print(f"Skipping trim_galore for {base_name} - already processed")
            return trimmed1, trimmed2
        
        cmd = f"trim_galore --quality 20 --paired --cores {cpu_count_safety_margin} --output_dir {output_dir} {fastq1} {fastq2}"
        subprocess.run(cmd, shell=True, check=True)
        self.save_checkpoint(base_name, "trim_galore")
        return trimmed1, trimmed2
    
    def run_bismark_alignment(self, trimmed1, trimmed2, reference_genome):
        """Align reads using Bismark"""
        output_dir = os.path.join(self.output_dir, "aligned")
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = os.path.basename(trimmed1).replace("_1_val_1.fq", "")
        bam_file = os.path.join(output_dir, f"{base_name}_bismark_bt2.bam")
        
        # Skip if file exists and stage is completed
        if os.path.exists(bam_file) and self.get_stage(base_name) == "bismark":
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


        cmd = f"bismark --parallel 8 --output_dir {output_dir} --genome {reference_genome} -1 {trimmed1} -2 {trimmed2}"
        subprocess.run(cmd, shell=True, check=True)
        self.save_checkpoint(base_name, "bismark")
        return bam_file
    
    def extract_methylation(self, bam_file, genes_of_interest, gene_coords):
        """Extract methylation data for specific genes"""
        output_dir = os.path.join(self.output_dir, "methylation")
        os.makedirs(output_dir, exist_ok=True)
        
        base_name = os.path.basename(bam_file).replace("_bismark_bt2.bam", "")
        methylation_file = os.path.join(output_dir, f"{base_name}_methylation.txt")
        
        # Skip if file exists and stage is completed
        if os.path.exists(methylation_file) and self.get_stage(base_name) == "methylation":
            print(f"Loading existing methylation data for {base_name}")
            return self.load_methylation_results(methylation_file)

        # Extract methylation data
        cmd = f"bismark_methylation_extractor --parallel {max([cpu_count_safety_margin // 3, 1])} --comprehensive --output {output_dir} {bam_file}"
        # NOTE: on --parallel:
        #> May also be --multicore <int>. Sets the number of cores to be used for the methylation extraction process.
        #> If system resources are plentiful this is a viable option to speed up the extraction process (we observed a
        #> near linear speed increase for up to 10 cores used). Please note that a typical process of extracting a BAM file
        #> and writing out '.gz' output streams will in fact use ~3 cores per value of --parallel <int>
        #> specified (1 for the methylation extractor itself, 1 for a Samtools stream, 1 for GZIP stream), so
        #> --parallel 10 is likely to use around 30 cores of system resources. This option has no bearing
        #> on the bismark2bedGraph or genome-wide cytosine report processes.
        subprocess.run(cmd, shell=True, check=True)
        
        # Analyze for each gene
        methylation_data = {}
        for gene in genes_of_interest:
            gene_bed = gene_coords[gene]
            cmd = f"bedtools intersect -a {methylation_file} -b {gene_bed}"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            methylation_data[gene] = self.parse_methylation_output(result.stdout)
        
        # Save results
        with open(methylation_file, 'w') as f:
            json.dump(methylation_data, f, indent=2)
        
        self.save_checkpoint(base_name, "methylation")
        return methylation_data
    
    @staticmethod
    def parse_methylation_output(methylation_text):
        """Parse methylation data and calculate statistics"""
        total_sites = 0
        methylated_sites = 0
        
        for line in methylation_text.split('\n'):
            if line:
                fields = line.split('\t')
                total_sites += 1
                if fields[3] == 'methylated':
                    methylated_sites += 1
        
        return {
            'total_sites': total_sites,
            'methylated_sites': methylated_sites,
            'methylation_percentage': (methylated_sites/total_sites*100) if total_sites > 0 else 0
        }
    
    @staticmethod
    def load_methylation_results(methylation_file):
        """Load previously calculated methylation results"""
        with open(methylation_file, 'r') as f:
            return json.load(f)

def analyze_methylation(sra_files, reference_genome, genes_of_interest, gene_coords):
    """Main analysis function"""
    analyzer = BisulfiteAnalyzer()
    results = {}
    
    for sra_file in sra_files:
        base_name = os.path.basename(sra_file)
        # print(f"\nProcessing {base_name}...")

        # # Convert SRA to FASTQ (if needed)
        # fastq1, fastq2 = analyzer.run_fastq_dump(sra_file)

        # print(f"Done with fastq conversion")
        # # Quality control and trimming (if needed)
        # trimmed1, trimmed2 = analyzer.run_trim_galore(fastq1, fastq2)
        # print(f"Done with Quality control and trimming")

        # # Alignment (if needed)
        # bam_file = analyzer.run_bismark_alignment(trimmed1, trimmed2, reference_genome)
        # TODO: fix that bam file is not returned correctly
        # TODO: make sure we have enough space for intermediate files
        bam_file = "methylation_analysis/aligned/SRR6228477_1_val_1_bismark_bt2_pe.bam"
        print(f"done with alignment")
        
        # Methylation analysis (if needed)
        results[base_name] = analyzer.extract_methylation(bam_file, genes_of_interest, gene_coords)
    
    return results

# Example usage
if __name__ == "__main__":
    genes_of_interest = ['NANOG', 'KLF17', 'DPPA4', 'NANOG_promoter', 'KLF17_promoter', 'DPPA4_promoter']
    
    # Define your file paths
    sra_files = ['/home/ubuntu/sra_data/sra/SRR6228477']
    reference_genome = '/home/ubuntu/sra_data/reference_genome/'
    gene_coords = {
        'NANOG': 'bedfiles/NANOG.bed',
        'KLF17': 'bedfiles/KLF17.bed',
        'DPPA4': 'bedfiles/DPPA4.bed',
        'NANOG_promoter': 'bedfiles/NANOG_promoter.bed',
        'KLF17_promoter': 'bedfiles/KLF17_promoter.bed',
        'DPPA4_promoter': 'bedfiles/DPPA4_promoter.bed',
    }

    assert len(genes_of_interest) == len(gene_coords)
    assert sorted(genes_of_interest) == sorted([x for x in gene_coords.keys()])
    
    # Run analysis
    results = analyze_methylation(sra_files, reference_genome, genes_of_interest, gene_coords)
    
    # Create summary DataFrame
    summary_df = pd.DataFrame.from_dict(results, orient='index')
    print(f"\nMethylation Analysis Summary (file:{results}):")
    print(summary_df)
