#!/usr/bin/env python

import sys

def process_bed(input_file):
    """Process BED file to filter out overlapping regions."""
    prev_line = None
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines
            if not line:
                continue
            parts = line.split('\t')
            # Skip invalid lines with less than 3 columns
            if len(parts) < 3:
                continue
            try:
                # Parse chromosome and coordinates
                chr_current = parts[0]
                start_current = int(parts[1])
                end_current = int(parts[2])
            except ValueError:
                # Skip lines with invalid coordinate values
                continue
            # First line processing
            if prev_line is None:
                print(line)
                prev_line = (chr_current, start_current, end_current)
            else:
                chr_prev, start_prev, end_prev = prev_line
                # Check if current region is completely within previous region
                if (chr_current == chr_prev and
                    start_current >= start_prev and
                    end_current <= end_prev):
                    # Skip this line as it's contained within previous region
                    continue
                else:
                    # Print line and update previous region
                    print(line)
                    prev_line = (chr_current, start_current, end_current)

if __name__ == '__main__':
    # Check command line arguments
    if len(sys.argv) != 2:
        print("Usage: python filter_bed.py <input.bed>")
        sys.exit(1)
    # Process the BED file
    process_bed(sys.argv[1])