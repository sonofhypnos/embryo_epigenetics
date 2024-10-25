#!/usr/bin/env python3

import pandas as pd
import requests
import os

def get_ensembl_coordinates(gene_name, species='human'):
    """
    Get gene coordinates from Ensembl REST API
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?"
    
    headers = {"Content-Type": "application/json"}
    r = requests.get(server + ext, headers=headers)
    
    if r.ok:
        decoded = r.json()
        return {
            'chromosome': decoded['seq_region_name'],
            'start': decoded['start'],
            'end': decoded['end'],
            'strand': decoded['strand'],
            'gene_name': gene_name
        }
    else:
        print(f"Failed to get coordinates for {gene_name}")
        return None

def create_bed_files(gene_list, output_dir='bed_files'):
    """
    Create BED files for a list of genes
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get coordinates for each gene and create BED files
    for gene in gene_list:
        print(f"Processing {gene}...")
        coords = get_ensembl_coordinates(gene)
        
        if coords:
            # Create BED file (BED format: chr start end name score strand)
            bed_path = os.path.join(output_dir, f"{gene}.bed")
            with open(bed_path, 'w') as f:
                # Note: BED format is 0-based, Ensembl is 1-based, so subtract 1 from start
                f.write(f"chr{coords['chromosome']}\t{coords['start']-1}\t{coords['end']}\t{gene}\t0\t{'+' if coords['strand'] == 1 else '-'}\n")
            
            # Add promoter region (2kb upstream)
            promoter_path = os.path.join(output_dir, f"{gene}_promoter.bed")
            with open(promoter_path, 'w') as f:
                if coords['strand'] == 1:  # Forward strand
                    promoter_start = max(0, coords['start'] - 2000)
                    f.write(f"chr{coords['chromosome']}\t{promoter_start-1}\t{coords['start']}\t{gene}_promoter\t0\t+\n")
                else:  # Reverse strand
                    promoter_start = coords['end']
                    f.write(f"chr{coords['chromosome']}\t{promoter_start}\t{promoter_start + 2000}\t{gene}_promoter\t0\t-\n")
            
            print(f"Created BED files for {gene}")
            print(f"Gene body: {bed_path}")
            print(f"Promoter region: {promoter_path}")
        else:
            print(f"Skipping {gene} - coordinates not found")

def extend_bed_regions(bed_file, upstream=2000, downstream=2000):
    """
    Extend regions in a BED file by specified amounts upstream and downstream
    """
    df = pd.read_csv(bed_file, sep='\t', header=None, 
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    for idx, row in df.iterrows():
        if row['strand'] == '+':
            df.at[idx, 'start'] = max(0, row['start'] - upstream)
            df.at[idx, 'end'] = row['end'] + downstream
        else:
            df.at[idx, 'start'] = max(0, row['start'] - downstream)
            df.at[idx, 'end'] = row['end'] + upstream
    
    output_file = bed_file.replace('.bed', f'_extended_{upstream}bp.bed')
    df.to_csv(output_file, sep='\t', header=False, index=False)
    return output_file

if __name__ == "__main__":
    # List of genes to process
    genes = ['NANOG', 'KLF17', 'DPPA4']
    
    # Create basic BED files
    create_bed_files(genes)
    
    # Optionally extend regions
    for gene in genes:
        bed_file = f"bed_files/{gene}.bed"
        if os.path.exists(bed_file):
            extended_file = extend_bed_regions(bed_file)
            print(f"Created extended region file: {extended_file}")
