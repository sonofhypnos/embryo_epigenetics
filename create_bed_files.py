#!/usr/bin/env python3

#!/usr/bin/env python3
import pandas as pd
import requests
import os
import json

class RegionMapping:
    def __init__(self):
        self.regions = {}  # Store original genomic coordinates
        self.index_map = {}  # Map custom index coordinates to genomic coordinates

    def add_region(self, chrom, start, end, name, custom_start=0):
        """
        Add a region to the mapping
        custom_start: The starting position in the custom index (usually 0)
        """
        region_length = end - start
        self.regions[name] = {
            'chrom': chrom,
            'orig_start': start,
            'orig_end': end,
            'custom_start': custom_start,
            'custom_end': custom_start + region_length
        }

        # Store mapping information
        self.index_map[name] = {
            'offset': start - custom_start,  # Offset to convert custom coords to genomic
            'chrom': chrom
        }

    def custom_to_genomic(self, name, custom_pos):
        """Convert a position in the custom index to genomic coordinates"""
        if name in self.index_map:
            return custom_pos + self.index_map[name]['offset']
        return None

    def save_mapping(self, output_file):
        """Save mapping information to a JSON file"""
        with open(output_file, 'w') as f:
            json.dump({
                'regions': self.regions,
                'index_map': self.index_map
            }, f, indent=2)

def get_ensembl_coordinates(gene_name, species='human'):
    """Get gene coordinates from Ensembl REST API"""
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

def create_bed_files_with_mapping(gene_list, output_dir='bed_files'):
    """Create BED files and coordinate mapping for a list of genes"""
    os.makedirs(output_dir, exist_ok=True)

    mapping = RegionMapping()
    created_files = []

    current_custom_start = 0  # Track position in custom index

    for gene in gene_list:
        print(f"Processing {gene}...")
        coords = get_ensembl_coordinates(gene)

        if coords:
            # Create standard gene body BED file
            bed_path = os.path.join(output_dir, f"{gene}.bed")
            with open(bed_path, 'w') as f:
                f.write(f"chr{coords['chromosome']}\t{coords['start']-1}\t{coords['end']}\t{gene}\t0\t{'+' if coords['strand'] == 1 else '-'}\n")
            created_files.append(bed_path)

            # Add to mapping
            mapping.add_region(
                f"chr{coords['chromosome']}",
                coords['start']-1,
                coords['end'],
                gene,
                current_custom_start
            )
            current_custom_start += (coords['end'] - (coords['start']-1))

            # Create promoter region BED file
            promoter_path = os.path.join(output_dir, f"{gene}_promoter.bed")
            with open(promoter_path, 'w') as f:
                if coords['strand'] == 1:
                    promoter_start = max(0, coords['start'] - 2000)
                    f.write(f"chr{coords['chromosome']}\t{promoter_start-1}\t{coords['start']}\t{gene}_promoter\t0\t+\n")
                    mapping.add_region(
                        f"chr{coords['chromosome']}",
                        promoter_start-1,
                        coords['start'],
                        f"{gene}_promoter",
                        current_custom_start
                    )
                else:
                    promoter_start = coords['end']
                    f.write(f"chr{coords['chromosome']}\t{promoter_start}\t{promoter_start + 2000}\t{gene}_promoter\t0\t-\n")
                    mapping.add_region(
                        f"chr{coords['chromosome']}",
                        promoter_start,
                        promoter_start + 2000,
                        f"{gene}_promoter",
                        current_custom_start
                    )
                current_custom_start += 2000
            created_files.append(promoter_path)

            print(f"Created BED files for {gene}")
        else:
            print(f"Skipping {gene} - coordinates not found")

    # Save the mapping
    mapping_file = os.path.join(output_dir, "coordinate_mapping.json")
    mapping.save_mapping(mapping_file)

    return created_files, mapping

def create_custom_index_regions(bed_file, buffer_size=500, output_dir='bed_files/custom_index'):
    """Create extended BED regions with buffer for custom index creation"""
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'name', 'score', 'strand'])

    # Add buffer to both ends
    df['start'] = df['start'].apply(lambda x: max(0, x - buffer_size))
    df['end'] = df['end'] + buffer_size

    # Create output filename
    base_name = os.path.basename(bed_file)
    output_file = os.path.join(output_dir, f"{os.path.splitext(base_name)[0]}_index.bed")

    # Save extended regions
    df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Created custom index region file: {output_file}")

    return output_file

def convert_methylation_coordinates(bedgraph_file, mapping_file, output_file):
    """Convert methylation coordinates from custom index to genomic coordinates"""
    # Load coordinate mapping
    with open(mapping_file, 'r') as f:
        mapping_data = json.load(f)

    # Read bedGraph file
    df = pd.read_csv(bedgraph_file, sep='\t', header=None,
                     names=['custom_chrom', 'custom_start', 'custom_end', 'methylation'])

    # Convert coordinates using mapping
    converted_rows = []
    for region_name, mapping_info in mapping_data['index_map'].items():
        offset = mapping_info['offset']
        chrom = mapping_info['chrom']

        # Filter rows for this region and convert coordinates
        region_rows = df[df['custom_chrom'] == region_name].copy()
        if not region_rows.empty:
            region_rows['custom_start'] = region_rows['custom_start'] + offset
            region_rows['custom_end'] = region_rows['custom_end'] + offset
            region_rows['custom_chrom'] = chrom
            converted_rows.append(region_rows)

    # Combine all converted coordinates
    if converted_rows:
        converted_df = pd.concat(converted_rows)
        converted_df.sort_values(['custom_chrom', 'custom_start'], inplace=True)
        converted_df.to_csv(output_file, sep='\t', header=False, index=False)
        return output_file
    return None

def process_all_regions(gene_list, buffer_size=500):
    """Process all genes and create both standard and custom index BED files"""
    # Create standard BED files with coordinate mapping
    standard_files, mapping = create_bed_files_with_mapping(gene_list)

    # Create custom index regions for each standard BED file
    custom_index_files = []
    for bed_file in standard_files:
        custom_file = create_custom_index_regions(bed_file, buffer_size)
        custom_index_files.append(custom_file)

    return standard_files, custom_index_files, mapping

if __name__ == "__main__":
    # List of genes to process
    genes = ['NANOG', 'KLF17', 'DPPA4']

    # Process all regions and create both standard and custom index files
    standard_files, custom_index_files, mapping = process_all_regions(genes, buffer_size=500)

    print("\nCreated standard BED files:")
    for f in standard_files:
        print(f"- {f}")

    print("\nCreated custom index BED files:")
    for f in custom_index_files:
        print(f"- {f}")
