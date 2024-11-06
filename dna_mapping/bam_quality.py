import pysam
from collections import defaultdict
import re
import sys

file, = sys.argv[1:]




def analyze_conversions(bam_file):
    """
    Analyze conversion patterns in bisulfite sequencing data.

    Parameters:
    bam_file (str): Path to BAM file
    """
    # Initialize counters
    conversions = defaultdict(int)
    total_mismatches = 0

    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    for read in bam.fetch():
        if read.has_tag("MD"):
            # Get the MD tag which shows the mismatches
            md = read.get_tag("MD")

            # Get the read sequence
            seq = read.query_sequence

            # Parse MD string to find positions and types of mismatches
            positions = re.findall(r"(\d+)([A-Z]|$)", md)
            current_pos = 0

            for length, ref_base in positions:
                current_pos += int(length)
                if ref_base and current_pos < len(seq):
                    read_base = seq[current_pos]
                    conversion = f"{ref_base}->{read_base}"
                    conversions[conversion] += 1
                    total_mismatches += 1
                current_pos += 1

    # Print results
    print("\nConversion Analysis Results:")
    print("-" * 40)
    print(f"Total mismatches analyzed: {total_mismatches}")
    print("\nConversion patterns:")
    for conversion, count in sorted(
        conversions.items(), key=lambda x: x[1], reverse=True
    ):
        percentage = (count / total_mismatches) * 100
        print(f"{conversion}: {count} ({percentage:.1f}%)")

    # Specifically look at C->T conversions
    ct_conversions = conversions.get("C->T", 0)
    ga_conversions = conversions.get("G->A", 0)  # Complement strand
    total_bisulfite_conversions = ct_conversions + ga_conversions

    print("\nBisulfite conversion analysis:")
    print(f"C->T conversions: {ct_conversions}")
    print(f"G->A conversions: {ga_conversions}")
    if total_mismatches > 0:
        print(
            f"Percentage of mismatches that are bisulfite conversions: {(total_bisulfite_conversions/total_mismatches)*100:.1f}%"
        )

    bam.close()


if __name__ == "__main__":
    # Replace with your BAM file path
    analyze_conversions(file)

