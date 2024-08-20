import os
import glob
import argparse
import numpy as np

#############
## Comments
##
## Currently inefficient for long sequences. Can be optimised by per-line scoring, instead of reading the whole sequence. See commented code
#############

# Function to parse arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Score input PWMs along the sequence in a fasta file')
    parser.add_argument('-i', '--path-to-fasta', nargs='+', required=True, help='Path(s) to input FASTA file with sequence(s) to score')
    parser.add_argument('-m', '--path-to-pwms', nargs='+', required=True, help='Path(s) to PWM(s) to score over input sequence(s)')
    parser.add_argument('-o', '--out-dir', default='.', help='Directory to write results (default: current directory)')
    parser.add_argument('-fs', '--suffix', type=str, default="", help='Desired suffix of the output bed file (default: none)')
    parser.add_argument('-s', '--pwm-score-lim', type=float, default=0.01, help='Lower limit of a pwm score to be recorded as a hit (default: 0.01)')
    args = parser.parse_args()
    expanded_pwms = []
    for pattern in args.path_to_pwms:
        expanded_pwms.extend(glob.glob(pattern))
    
    if not expanded_pwms:
        raise FileNotFoundError("No PWM files found matching the provided patterns.")
    
    args.path_to_pwms = expanded_pwms
    return args

# Read PWM from a file
def read_pwm(pwm_path):
    with open(pwm_path) as pwm_file:
        pwm = [line.strip().split('\t')[1:] for line in pwm_file.readlines()]
        bp_order = pwm[0]
        pwm = np.array(pwm[1:], dtype=float)
    return pwm, bp_order

# Read fasta file into a dictionary
def read_fasta_as_dict(fasta_path):
    """Read FASTA file into a dictionary."""
    sequences = {}
    with open(fasta_path) as fasta_file:
        seq_name = None
        sequence = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if seq_name:
                    sequences[seq_name] = ''.join(sequence).replace('T', 'U')
                seq_name = line[1:].replace(":", "\t").replace("-", "\t")
                sequence = []
            else:
                sequence.append(line)
        if seq_name:
            sequences[seq_name] = ''.join(sequence).replace('T', 'U')
    return sequences

# Score a given fasta sequence with PWM
def score_sequence(seq, pwm, bp_order, score_lim):
    mot_l = pwm.shape[0]
    seq_length = len(seq)
    n_iter = seq_length - mot_l + 1
    scores = []
    for i in range(n_iter):
        subseq = seq[i:i + mot_l]
        prod = 1.0
        min_prod = 1.0
        max_prod = 1.0
        
        for j, base in enumerate(subseq):
            try:
                idx = bp_order.index(base)
            except ValueError:
                prod = 0.0
                min_prod = 0.0
                max_prod = 0.0
                break
            
            probs = pwm[j]
            prob = probs[idx]
            max_prob = np.max(probs)
            min_prob = np.min(probs)
            
            prod *= prob
            min_prod *= min_prob
            max_prod *= max_prob

        if max_prod != min_prod:
            mss = (prod - min_prod) / (max_prod - min_prod)
            if mss > score_lim:
                scores.append((i, i + mot_l, mss))
    
    return scores

# Main
def main():
    args = parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    if args.suffix != "":
        suffix = "_" + args.suffix
    else:
        suffix = args.suffix
        
    for fasta_path in args.path_to_fasta:
            fasta_name = os.path.basename(fasta_path).replace('.fasta', '').replace('.fa', '')
            output_file = os.path.join(args.out_dir, f"{fasta_name}{suffix}_scored.bed")

            sequences = read_fasta_as_dict(fasta_path)

            with open(output_file, 'w') as out_f:
                for seq_name, sequence in sequences.items():
                    for pwm_path in args.path_to_pwms:
                        pwm, bp_order = read_pwm(pwm_path)
                        pwm_name = os.path.basename(pwm_path).replace('.txt', '')
                        scores = score_sequence(sequence, pwm, bp_order, args.pwm_score_lim)
                        
                        for start, end, mss in scores:
                            out_f.write(f"{seq_name}\t{start}\t{end}\t{pwm_name}\t{mss}\t.\n")

if __name__ == '__main__':
    main()


# def main():
#     args = parse_args()

#     if not os.path.exists(args.out_dir):
#         os.makedirs(args.out_dir)

#     for fasta_path in args.path_to_fasta:
#         fasta_name = os.path.basename(fasta_path).replace('.fasta', '').replace('.fa', '')
#         output_file = os.path.join(args.out_dir, f"{fasta_name}_scored.bed")

#         with open(output_file, 'w') as out_f:
#             with open(fasta_path) as fasta_file:
#                 for line in fasta_file:
#                     line = line.strip()
#                     if line.startswith(">"):
#                         seq_name = line[1:].replace(":", "\t").replace("-", "\t")
#                         seq = next(fasta_file).strip().replace('T', 'U')
                        
#                         for pwm_path in args.path_to_pwms:
#                             pwm, bp_order = read_pwm(pwm_path)
#                             pwm_name = os.path.basename(pwm_path).replace('.txt', '')
#                             scores = score_sequence(seq, pwm, bp_order)
                            
#                             for start, end, mss in scores:
#                                 out_f.write(f"{seq_name}\t{start}\t{end}\t{pwm_name}\t{mss}\t.\n")

# if __name__ == "__main__":
#     main()