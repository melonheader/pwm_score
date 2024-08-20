# pwm_score

## Foreword
This repository contains a small toolkit designed to swiftly score arbitrary position-weight matrices (PWMs) against a DNA sequence in FASTA format. Included in the repository is a collection of PWMs for known RNA-binding proteins (RBPs).

## Usage

### Extracting the PWM Collection
If you need to extract the PWM collection, use the following command:
```bash
tar -xzf data.tar.gz
```
### Main tool
The package contains one main function for the scoring - `digest.py`:
```bash
python pwm_score/digest.py --help
usage: digest.py [-h] -i PATH_TO_FASTA [PATH_TO_FASTA ...] -m PATH_TO_PWMS [PATH_TO_PWMS ...] [-o OUT_DIR]
                 [-fs SUFFIX] [-s PWM_SCORE_LIM]

Score input PWM along the sequence

options:
  -h, --help            show this help message and exit
  -i PATH_TO_FASTA [PATH_TO_FASTA ...], --path-to-fasta PATH_TO_FASTA [PATH_TO_FASTA ...]
                        Path(s) to input FASTA file with sequence(s) to score
  -m PATH_TO_PWMS [PATH_TO_PWMS ...], --path-to-pwms PATH_TO_PWMS [PATH_TO_PWMS ...]
                        Path(s) to PWM(s) to score over input sequence(s)
  -o OUT_DIR, --out-dir OUT_DIR
                        Directory to write results (default: current directory)
  -fs SUFFIX, --suffix SUFFIX
                        Desired suffix of the output bed file (default: none)
  -s PWM_SCORE_LIM, --pwm-score-lim PWM_SCORE_LIM
                        Lower limit of a pwm score to be recorded as a hit (default: 0.01)
```
### Example run
The standard use case implies a fasta file with an arbitrary number of entries and a selection of PWMs one want to score. 
For example, the `test` directory contains a fasta seuqence `fasta_test.fa` with two entries, first of which carries a perfect match to the motif from `M112_0.6.txt` PWM.
Running the main pacakge's function `digest.py` over the `fasta_test.fa` with all PWMs in the `data/pwm/` folder as follows
```bash
python digest.py -i test/fasta_test.fa -m test/pwms/*pwm -o test
```
produces an output test/fasta_test_scored.bed in the bed format (0-based) where the names of fasta entries are recorded in the first column,
hit coordinates are in the columns 2 and 3, and the name of the PWM file corresponding to the hit and the hit's score in the columns 4 and 5.

Finally, let's validate the scoring by checking for the score of an artificially inserted perfect match: 
```bash
head -3 test/fasta_test_scored.bed
ayy_1	10	16	motif1.pwm	0.09777396784654607	.
ayy_1	28	34	motif2.pwm	**1.0**	.
ayy_1	29	35	motif2.pwm	0.020579816903596112	.
```

### Visualisation

The scored bed files can be further digested into wigs to be visualised in a genomic browser of your choice ny the `digest_to_wig.py` function:
```bash
python pwm_score/digest_to_wig.py --help
usage: digest_to_wig.py [-h] -b BED_FILE [-o OUTPUT_DIR] [-sl SPAN_LENGTH] -so START_OFFSET

Convert BED file to WIG files split by PWM names.

options:
  -h, --help            show this help message and exit
  -b BED_FILE, --bed_file BED_FILE
                        Path to the input BED file.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to write WIG files.
  -sl SPAN_LENGTH, --span-length SPAN_LENGTH
                        Length of the PWM record to adjust the span.
  -so START_OFFSET, --start_offset START_OFFSET
                        Starting genomic coordinate (chrN:position).
```
The `-sl` parameter controls the span of the motif hit in bases to be recorded in the wig track. In our case, it is set to 7 as all motifs are 7 bases long.
The `-so` parameter is needed to adjust the relative coordinates of the bed file with genomic coordinates.

Running the function over the scored test bed file
```bash
python pwm_score/digest_to_wig.py -b test/fasta_test_scored.bed -o test -sl 7 -so chr2:1337
```
produces a separate wig file per motif that has at least one hit.