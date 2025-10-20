#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Filter file lines',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('-d', '--down', type=float)
    parser.add_argument('-u', '--upper', type=float)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate arguments
    if args.down >= args.upper:
        print("Error: Lower bound (-d) must be less than upper bound (-u)")
        sys.exit(1)
    
    try:
        with open(args.input_file, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        print("File read error: {}".format(e))
        sys.exit(1)
    
    # Output first line
    if len(lines) > 0:
        print(lines[0].rstrip('\n'))
    
    # Exit if only one line exists
    if len(lines) <= 1:
        return
    
    # Store previous valid line for comparison
    prev_line = lines[0].strip().split('\t')
    
    # Process from second line onward
    for i in range(1, len(lines)):
        current_line = lines[i].strip().split('\t')
        
        # Ensure line has at least 6 columns
        if len(current_line) < 6 or len(prev_line) < 6:
            continue
        
        try:
            # Extract values from column 2 (index 1) and column 6 (index 5)
            curr_col2 = int(current_line[1])
            prev_col2 = int(prev_line[1])
            curr_col6 = int(current_line[5])
            prev_col6 = int(prev_line[5])
            
            number1 = curr_col2 - prev_col2 + 1
            number2 = curr_col6 - prev_col6 + 1
            
            if number2 == 0:
                continue
            
            number3 = float(number1) - float(number2)
            
            # Check if condition is met (using command line arguments)
            if args.down <= number3 <= args.upper:
                # Output current line
                print(lines[i].rstrip('\n'))
                # Update prev_line for next comparison
                prev_line = current_line
            # Skip line if condition not met (prev_line not updated)
            
        except (ValueError, IndexError) as e:
            # Skip lines with invalid data
            continue

if __name__ == "__main__":
    main()