import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Convert BED file to WIG files split by PWM names.')
    parser.add_argument('-b', '--bed_file', required=True, help='Path to the input BED file.')
    parser.add_argument('-o', '--output_dir', default='.', help='Directory to write WIG files.')
    parser.add_argument('-sl', '--span-length', default=1, help='Length of the PWM record to adjust the span.')
    parser.add_argument('-so', '--start_offset', required=True, help='Starting genomic coordinate (chrN:position).')
    
    args = parser.parse_args()
    return args

def parse_bed_file(bed_path):
    """Parse BED file and group scores by PWM name."""
    pwm_dict = {}
    with open(bed_path) as bed_file:
        for line in bed_file:
            if line.strip():  # Skip empty lines
                columns = line.strip().split()
                if len(columns) < 6:
                    continue  # Skip invalid lines
                chrom, start, end, pwm_name, score, _ = columns
                start, end, score = int(start), int(end), float(score)
                
                if pwm_name not in pwm_dict:
                    pwm_dict[pwm_name] = []
                pwm_dict[pwm_name].append((start, end, score))
    return pwm_dict

def write_wig_file(pwm_name, scores, user_chrom, start_offset, span_length, output_dir):
    """Write scores to WIG file."""
    output_path = os.path.join(output_dir, f"{pwm_name}.wig")
    with open(output_path, 'w') as wig_file:
        wig_file.write(f"track type=wiggle_0 name=\"{pwm_name}\"\n")
        wig_file.write(f"variableStep chrom={user_chrom} span={span_length}\n")
        for start, end, score in scores:
            adjusted_start = start_offset + start
            wig_file.write(f"{adjusted_start}\t{score}\n")

def main():
    args = parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Parse start offset
    user_chrom, start_pos = args.start_offset.split(':')
    start_offset = int(start_pos)

    # Parse BED file
    pwm_dict = parse_bed_file(args.bed_file)

    # Write WIG files
    
    for pwm_name, scores in pwm_dict.items():
        write_wig_file(pwm_name, scores, user_chrom, start_offset, args.span_length, args.output_dir)

if __name__ == '__main__':
    main()