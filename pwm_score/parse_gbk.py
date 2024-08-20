import argparse
from Bio import SeqIO
import re

def parse_args():
    parser = argparse.ArgumentParser(description="Convert a GenBank (.gbk) file to FASTA and BED format.")
    parser.add_argument("-gbk", "--gbk-file", type=str, help="Path to the input GenBank (.gbk) file")
    parser.add_argument("-fo", "--fasta-output", type=str, help="Path to the output FASTA file")
    parser.add_argument("-bo", "--bed-output", type=str, help="Path to the output BED file")
    args = parser.parse_args()
    return args

def parse_location(location):
    """
    Parses the GenBank location field to handle 'join' and 'complement' formats.
    Returns a list of (start, end) tuples.
    """
    ranges = []

    # Remove any complement() if present
    location = re.sub(r'complement\(|\)', '', location)
    # Handle 'join' and split into individual ranges
    if 'join' in location:
        location = location.replace('join{', '').replace('}', '')
        # Extract the range segments
        segments = re.split(r',\s*', location)
        for segment in segments:
            segment = segment.strip().replace('[', '').replace('](+', '')
            if ':' in segment:
                try:
                    start, end = map(int, segment.split(':'))
                    ranges.append((start, end))
                except ValueError:
                    print(f"Warning: Unable to parse segment '{segment}'")
    else:
        # Handle non-join ranges
        segment = location.strip()
        if '..' in segment:
            try:
                start, end = map(int, segment.split('..'))
                ranges.append((start, end))
            except ValueError:
                print(f"Warning: Unable to parse segment '{segment}'")
    
    if not ranges:
        print(f"Warning: No valid ranges found for location '{location}'")
    
    return ranges

def gbk_to_fasta_and_bed(gbk_file, fasta_output, bed_output):
    with open(gbk_file, "r") as gb_file:
        record = SeqIO.read(gb_file, "genbank")
        
        # Write the sequence to a FASTA file
        with open(fasta_output, "w") as fasta_file:
            SeqIO.write(record, fasta_file, "fasta")
        
        # Write the features to a BED file
        with open(bed_output, "w") as bed_file:
            for feature in record.features:
                # Check for the /label qualifier first
                if "label" in feature.qualifiers:
                    feature_name = feature.qualifiers["label"][0]
                elif "gene" in feature.qualifiers:
                    feature_name = feature.qualifiers["gene"][0]
                elif "locus_tag" in feature.qualifiers:
                    feature_name = feature.qualifiers["locus_tag"][0]
                else:
                    feature_name = "unknown"
                
                # Handle different types of features
                if feature.type in ["mRNA", "CDS"]:
                    # Parse location field
                    locations = parse_location(str(feature.location))
                    
                    # Write each segment to the BED file
                    for start, end in locations:
                        strand = "+" if feature.location.strand == 1 else "-"
                        bed_file.write(f"{record.id}\t{start}\t{end}\t{feature_name}\t0\t{strand}\n")
                else:
                    # Write other features normally
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"
                    bed_file.write(f"{record.id}\t{start}\t{end}\t{feature_name}\t0\t{strand}\n")
def main():
    args = parse_args()
    gbk_to_fasta_and_bed(args.gbk_file, args.fasta_output, args.bed_output)
    print(f"FASTA file written to: {args.fasta_output}")
    print(f"BED file written to: {args.bed_output}")

if __name__ == "__main__":
    main()