#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Number Sequence Generator
This script generates a sequence of 9 numbers where each number multiplied by 2 equals the next number,
and the 9th number multiplied by 2 equals the input value.
"""

import sys
import argparse


def generate_sequence(target_value):
    """
    Generate a sequence of 9 numbers where each number doubled equals the next,
    and the last number doubled equals the target value.
    
    Args:
        target_value: The final target number (9th number * 2)
    
    Returns:
        A list of 9 integers forming the sequence (rounded if necessary)
    """
    # The 9th number should be target_value / 2
    # The 8th number should be target_value / 4
    # The 7th number should be target_value / 8
    # ...
    # The 1st number should be target_value / 512 (2^9)
    
    sequence = []
    divisor = 2 ** 9  # 512
    
    for i in range(9):
        number = target_value / float(divisor)
        # Round to nearest integer
        rounded_number = int(round(number))
        sequence.append(rounded_number)
        divisor = divisor / 2
    
    return sequence


def main():
    """Main function to parse arguments and generate the sequence"""
    parser = argparse.ArgumentParser(
        description='Generate a sequence of 9 numbers where each number multiplied by 2 equals the next number, '
                    'and the 9th number multiplied by 2 equals the input value.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -n 3000000
  %(prog)s --number 1000000
        """
    )
    
    parser.add_argument(
        '-n', '--number',
        type=float,
        required=True,
        help='The target number (the value that the 9th number multiplied by 2 should equal)'
    )
    
    args = parser.parse_args()
    
    target = args.number
    
    # Generate the sequence
    sequence = generate_sequence(target)
    
    # Print the sequence (one number per line)
    for num in sequence:
        print num
    
    # Print the final target value (rounded to integer)
    print int(round(target))


if __name__ == '__main__':
    main()
