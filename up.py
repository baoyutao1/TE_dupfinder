import argparse

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="BED file proximity filter")
    parser.add_argument("input", help="Input BED file path")
    parser.add_argument("-c", "--coordinate", type=int, required=True,
                        help="Target genomic coordinate to check proximity (required)")
    # Add distance parameter with default value 500
    parser.add_argument("-d", "--distance", type=int, default=500,
                        help="Proximity distance threshold (default: 500)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    with open(args.input, 'r') as f:
        lines = [line.strip().split() for line in f]

    if not lines:
        return

    # Process first line
    first_line = lines[0]
    if len(first_line) < 4:
        return

    try:
        chr1, start1, end1 = first_line[0], int(first_line[1]), int(first_line[2])
    except (IndexError, ValueError):
        return

    # Coordinate proximity check using distance parameter
    if abs(args.coordinate - end1) > args.distance:
        return

    output = [first_line]
    prev_start, prev_end = start1, end1

    for line in lines[1:]:
        if len(line) < 4:
            continue

        if line[0] != chr1:  # Stop if chromosome changes
            break

        try:
            curr_start = int(line[1])
            curr_end = int(line[2])
        except (IndexError, ValueError):
            continue

        # Calculate overlap or gap
        overlap = curr_start < prev_end and curr_end > prev_start
        gap = max(0, curr_start - prev_end) if curr_start >= prev_end else max(0, prev_start - curr_end)
        
        # Use distance threshold for proximity check
        if overlap or gap <= args.distance:
            output.append(line)
            prev_start, prev_end = curr_start, curr_end
        else:
            break

    # Print results in Python 2 compatible way
    for line in output:
        print '\t'.join(map(str, line))

if __name__ == "__main__":
    main()