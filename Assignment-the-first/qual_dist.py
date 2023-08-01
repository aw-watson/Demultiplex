#!/usr/bin/env python

#Generate a per base distribution of quality scores for read1, read2, index1, and index2. 
#Average the quality scores at each position for all reads and generate a per nucleotide mean distribution.

import argparse
from matplotlib import pyplot as plt
import gzip


def convert_phred(letter: str) -> int:
    '''Converts a single character (Phred+33 encoding) into a numerical phred score'''
    return ord(letter)-33

def get_args():
    parser = argparse.ArgumentParser(description = "Calculate per-base position average quality scores within a gzipped FASTQ file")
    parser.add_argument("-f", "--filename", help="FASTQ file with ", required = True)
    parser.add_argument("-o", "--output-file", help="Name of output file to save created plot to. Will be appended with .png", required = True)
    parser.add_argument("-l", "--label", help="Description of input file. Used to label plots.")
    return parser.parse_args()

#main script

#read in command-line arguments, assign to variables
args = get_args()
fq: str = args.filename
lbl: str = args.label
outfile: str = args.output_file
#variable setup
qlines:int = 0 #for ease of calculation later
means: list = []

with gzip.open(fq, 'rt') as f:
    f.readline() #skip header
    f.readline() #skip sequence
    f.readline() #skip + line
    sequence: str = f.readline().strip() #first quality line
    qlines += 1
    for score in sequence: #set up list
        means.append(convert_phred(score))
    while True:
        f.readline()
        f.readline()
        f.readline()
        sequence = f.readline().strip()
        if not sequence:
            break
        qlines += 1
        for i, score in enumerate(sequence):
            means[i] += convert_phred(score)

for i, sum in enumerate(means):
    means[i] = sum/qlines

#plot
plt.rcParams["figure.figsize"] = (30,15)
plt.rcParams["font.size"] = 20
plt.title(f"Mean Phred Score by Base Position for {lbl}")
plt.xlabel("Base Position")
plt.ylabel("Average Phred Score")
plt.yscale('linear')
plt.bar(range(len(means)), means)
plt.savefig(f"{outfile}.png")


    