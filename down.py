import sys
import argparse

def process_bed(input_file, distance=500):
    """
    Process BED file to filter regions based on proximity.
    
    Args:
        input_file (str): Path to input BED file
        distance (int): Maximum allowed gap between regions (default: 500)
    
    Returns:
        list: Filtered BED lines
    """
    with open(input_file, 'r') as f:
        # Read non-empty lines
        lines = [line.strip() for line in f if line.strip()]

    if not lines:
        return []

    # Process first line
    first_line = lines[0]
    parts = first_line.split('\t')
    # Skip if less than 4 columns
    if len(parts) < 4:
        return []
    
    try:
        chrom = parts[0]
        start = int(parts[1])
    except ValueError:
        return []
    
    # First proximity check using distance parameter
    if start - 1 > distance:
        return []
    
    # Prepare output
    output = [first_line]
    prev_end = int(parts[2])
    
    # Process subsequent lines
    for line in lines[1:]:
        parts = line.split('\t')
        # Skip lines with less than 4 columns
        if len(parts) < 4:
            continue
        
        try:
            curr_start = int(parts[1])
            curr_end = int(parts[2])
        except ValueError:
            continue
        
        # Check proximity using distance parameter
        if curr_start <= prev_end + distance:
            output.append(line)
            prev_end = curr_end
        else:
            break
    
    return output

def main():
    """Main function to handle command line arguments and output."""
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Filter BED regions based on proximity')
    parser.add_argument('input', help='Input BED file path')
    parser.add_argument('-d', '--distance', type=int, default=500,
                        help='Maximum allowed gap between regions (default: 500)')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process BED file
    result = process_bed(args.input, args.distance)
    
    # Output results
    for line in result:
        print line

if __name__ == '__main__':
    main()