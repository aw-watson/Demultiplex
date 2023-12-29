#!/usr/bin/env python
#Author: André Watson, August 2023

import typing #for annotations
import bioinfo
import os
import argparse
import re
import gzip
import numpy as np

def amend_headers(R1: list, R2: list, R3: list, R4: list) -> None:
    '''Takes four properly formatted FASTQ records and adds the sequence line of R2 and R3 
    to the headers of both R1 and R4.'''
    R1[0] += ":" + R2[1] + "-" + R3[1]
    R4[0] += ":" + R2[1] + "-" + R3[1]
    return

def get_args():
    parser = argparse.ArgumentParser(description = "Demultiplex Illumina sequencing data")
    parser.add_argument("-i", "--input-directory", help="Directory with 4 gzipped FASTQ files containing Illumina sequencing data. \
                        Expects files to be identically named except for R[1|2|3|4] in the filename. \
                        R1 is the first read, R2 is the first index read, R3 is the second index read, R4 is the second read", required = True)
    parser.add_argument("-o", "--output-directory", help="Directory to save demultiplexed output files to.", required = True)
    parser.add_argument("-x", "--index-list", help="Text file with list of indices to demultiplex. \
                        Expects tab-separated data, with the first field as the index code, and the second field as the sequence", required = True)
    parser.add_argument("-q", "--quality-cutoff", help = "If an index has an average Phred score below this value, classify it as unknown", type = int, required = True)
    return parser.parse_args()

#script body

#assign arguments to variables--partly for convenience, partly for annotations
args = get_args()
in_list : list = os.listdir(args.input_directory)
out_dir: str = args.output_directory
idx_table: str = args.index_list
qual_cutoff: int = args.quality_cutoff

#extract naming convention from input files: to be used in naming output files
prefix: str = ""
suffix: str = ""
for file in in_list: #have to check within the directory for our input files
    m = re.fullmatch("(.*)(R[1234])(.*)\\.fastq.gz", file)
    if m:
        prefix = m.group(1)
        suffix = m.group(3)
    else:
        in_list.remove(file)

#put known indices in a dictionary and track occurrences
match_ctr: dict[str,int] = {}
with open(idx_table, 'rt') as idx_f:
    for line in idx_f:
        split_line = line.strip().split('\t')
        #map index code to occurrences, initialized at 0
        match_ctr[split_line[1]] = 0

#initialize mismatch dictionary and unknown counter
mism_ctr: dict[str, int] ={}
bad_idx_ctr: int = 0

#initialize lists to hold a FASTQ record at a time
r1: list = ["","","",""]
r2: list = ["","","",""]
r3: list = ["","","",""]
r4: list = ["","","",""]

#open input files, storing file handles in dictionary
input_files: dict[str, typing.TextIO] = {}
for filename in in_list:
    input_files[filename] = gzip.open(f"{args.input_directory}/{filename}", 'rt')

#open output files, storing file handles in dictionary
output_files: dict[str, typing.TextIO] = {}
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
#open output files for R1 and R2 with matching indices
for sequence in match_ctr.keys():
    output_files[f"R1_{sequence}"] = open(f"{out_dir}/{prefix}R1{suffix}_{sequence}.fastq", 'wt')
    output_files[f"R2_{sequence}"] = open(f"{out_dir}/{prefix}R2{suffix}_{sequence}.fastq", 'wt')
#open output files for R1 and R2 for mismatched and unknown index pairs
output_files["R1_mismatched"] = open(f"{out_dir}/{prefix}R1{suffix}_mismatched.fastq",'wt')
output_files["R2_mismatched"] = open(f"{out_dir}/{prefix}R2{suffix}_mismatched.fastq",'wt')
output_files["R1_unknown"] = open(f"{out_dir}/{prefix}R1{suffix}_unknown.fastq",'wt')
output_files["R2_unknown"] = open(f"{out_dir}/{prefix}R2{suffix}_unknown.fastq",'wt')

#main loop
while True:
    #as os.listdir returns in arbitrary order, we can't loop over input_list to get our input filenames
    #instead, reconstruct each filename to take a record from, and loop over the four lines of a FASTQ record
    for i in range(4):
        r1[i] = input_files[f"{prefix}R1{suffix}.fastq.gz"].readline().strip()
        r2[i] = input_files[f"{prefix}R2{suffix}.fastq.gz"].readline().strip()
        r3[i] = input_files[f"{prefix}R3{suffix}.fastq.gz"].readline().strip()
        r4[i] = input_files[f"{prefix}R4{suffix}.fastq.gz"].readline().strip()
    #if we've consumed all FASTQ records, and readline only returns empty strings, we're done
    if not r1[0]:
        break
    r3[1] = bioinfo.reverse_complement(r3[1])  # type: ignore
    amend_headers(r1,r2,r3,r4)

    #if indices are unknown, or if they're below a quality cutoff, put into the unknown files
    if (r2[1] not in match_ctr or 
            r3[1] not in match_ctr or 
            bioinfo.qual_score(r2[3]) < qual_cutoff or  # type: ignore
            bioinfo.qual_score(r3[3]) < qual_cutoff):  # type: ignore
        #increment counter
        bad_idx_ctr += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files["R1_unknown"].writelines(r1)
        output_files["R2_unknown"].writelines(r4)
    #else if we have valid but mismatched indices
    elif r2[1] != r3[1]:
        #increment counter
        idx_pair = r2[1] + "\t" + r3[1]
        if idx_pair not in mism_ctr:
            mism_ctr[idx_pair] = 1
        else:
            mism_ctr[idx_pair] += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files["R1_mismatched"].writelines(r1)
        output_files["R2_mismatched"].writelines(r4)
    #else if we have matching indices
    else: 
        #increment counter
        idx: str = r2[1]
        match_ctr[idx] += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files[f"R1_{idx}"].writelines(r1)
        output_files[f"R2_{idx}"].writelines(r4)
#close input files
for file in input_files.values():
    file.close()
#close output files
for file in output_files.values():
    file.close()

#summary information
mismatch_sum: int = 0
total_matched: int = sum(match_ctr.values())

#mismatched index pairs, output as a table
with open(f"{out_dir}/mismatched_counts.tsv", 'wt') as mmfh:
    mismatch_holder = np.zeros((24,24))
    codes = list(match_ctr.keys()) #ordered list
    for mm_pair in mism_ctr.keys():
        mismatch_sum += mism_ctr[mm_pair]
        mismatch = mm_pair.split("\t")
        mismatch_holder[codes.index(mismatch[0])][codes.index(mismatch[1])] = mism_ctr[mm_pair]
    mmfh.write("i,j\t" + "\t".join(codes) + "\n") #header row (column labels)
    for i, code in enumerate(codes):
        mmfh.write(code + '\t') #row label
        for v in mismatch_holder[i]:
            mmfh.write(str(v)+'\t') 
        mmfh.write('\n')

with open(f"{out_dir}/matched_counts.tsv", 'wt') as mfh:
    mfh.write("Index\tOccurrences\tPercent of Successfully Demultiplexed Reads\n")
    for idx_pair in match_ctr.keys():
        mfh.write(f"{idx_pair}\t{match_ctr[idx_pair]}\t{round(match_ctr[idx_pair]/total_matched*100, 2)}\n")
        if match_ctr[idx_pair] == 0:
            os.remove(f"{out_dir}/{prefix}R1{suffix}_{idx_pair}.fastq")
            os.remove(f"{out_dir}/{prefix}R2{suffix}_{idx_pair}.fastq")

print(f"Number of unknown or low-quality index pairs:\t{bad_idx_ctr}\n\
    Total number of mismatched index pairs:\t{mismatch_sum}\n\
    Total number of matching index pairs:\t{total_matched}")


