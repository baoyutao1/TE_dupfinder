import sys

def merge_bed_regions(input_file):
    """Merge all regions in a BED file into a single region."""
    chromosomes = []
    starts = []
    ends = []
    
    with open(input_file, 'r') as f:
        for line_number, line in enumerate(f, 1):
            line = line.strip()
            # Skip empty lines
            if not line:
                continue
            
            parts = line.split()
            # Check for minimum required columns (chrom, start, end)
            if len(parts) < 3:
                sys.stderr.write("Error: Line %d has insufficient columns (minimum 3 required)\n" % line_number)
                return None
            
            try:
                # Parse chromosome and coordinates
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                sys.stderr.write("Error: Line %d contains non-numeric coordinates\n" % line_number)
                return None
            
            # Store values for processing
            chromosomes.append(chrom)
            starts.append(start)
            ends.append(end)
    
    # Check if any valid data was found
    if not chromosomes:
        sys.stderr.write("Error: File is empty or contains no valid data\n")
        return None
    
    # Verify all regions are on the same chromosome
    if len(set(chromosomes)) != 1:
        sys.stderr.write("Error: Multiple chromosomes detected\n")
        return None
    
    # Return merged region (chromosome, min start, max end)
    return (chromosomes[0], min(starts), max(ends))

def main():
    """Main function to handle command line arguments and output."""
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python merge_bed.py input.bed\n")
        sys.stderr.write("Example: python merge_bed.py input.bed > merged.bed\n")
        sys.exit(1)
    
    # Process the BED file
    result = merge_bed_regions(sys.argv[1])
    if result:
        # Output merged region in BED format
        print "%s\t%d\t%d" % (result[0], result[1], result[2])

if __name__ == "__main__":
    main()