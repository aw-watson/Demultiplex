#!/usr/bin/env python
#Author: André Watson, August 2023

import typing #for annotations
import bioinfo
import os
import argparse
import re
import gzip

def amend_headers(R1: list, R2: list, R3: list, R4: list) -> None:
    '''Takes four properly formatted FASTQ records and adds the sequence line of R2 and R3 
    to the headers of both R1 and R4.'''
    R1[0] += ":" + R2[1] + "-" + R3[1]
    R4[0] += ":" + R2[1] + "-" + R3[1]
    return

def construct_mismatch(idx1: str, idx2:str) -> str:
    '''Helper function. Given two strings, returns them in alphabetical order, separated by a hyphen.
    Therefore, returns one string per combination of inputs, regardless of order'''
    if idx1 < idx2:
        return idx1 + "-" + idx2
    else:
        return idx2 + "-" + idx1


def get_args():
    parser = argparse.ArgumentParser(description = "Demultiplex Illumina sequencing data")
    parser.add_argument("-i", "--input-directory", help="Directory with 4 gzipped FASTQ files containing Illumina sequencing data. \
                        Expects files to be identically named except for R[1|2|3|4] in the filename. \
                        R1 is the first read, R2 is the first index read, R3 is the second index read, R4 is the second read", required = True)
    parser.add_argument("-o", "--output-directory", help="Directory to save demultiplexed output files to.", required = True)
    parser.add_argument("-x", "--index-list", help="Text file with list of indices to demultiplex. \
                        Expects tab-separated data, with the first field as the index code, and the second field as the sequence", required = True)
    return parser.parse_args()

#script body

#TODO: input validation/safeguards

#assign arguments to variables--partly for convenience, partly for annotations
args = get_args()
in_list : list = os.listdir(args.input_directory)
out_dir: str = args.output_directory
idx_table: str = args.index_list

#extract naming convention from input files: to be used in naming output files
prefix: str = ""
suffix: str = ""
m = re.fullmatch("(.*)(R[1234])(.*)\\.fastq.gz", in_list[0])
if m:
    prefix = m.group(1)
    suffix = m.group(3)
else:
    raise ValueError("Expected only gzipped fastq files as input, with R1, R2, R3, or R4 in the filename.")

#put known indices in a dictionary and track occurrences
idx_dict: dict = {}
match_ctr: dict = {}
with open(idx_table, 'rt') as idx_f:
    for line in idx_f:
        split_line = line.strip().split('\t')
        #map sequence to index code
        idx_dict[split_line[1]] = split_line[0]
        #map index code to occurrences, initialized at 0
        match_ctr[split_line[0]] = 0

#initialize mismatch dictionary and unknown counter
mism_ctr: dict ={}
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
#R1 and R2 with matching indices
for code in idx_dict.values():
    output_files[f"{prefix}R1{suffix}_{code}.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R1{suffix}_{code}.fastq.gz", 'wt')
    output_files[f"{prefix}R2{suffix}_{code}.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R2{suffix}_{code}.fastq.gz", 'wt')
#R1 and R2 for mismatched and unknown
output_files[f"{prefix}R1{suffix}_mismatched.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R1{suffix}_mismatched.fastq.gz",'wt')
output_files[f"{prefix}R2{suffix}_mismatched.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R2{suffix}_mismatched.fastq.gz",'wt')
output_files[f"{prefix}R1{suffix}_unknown.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R1{suffix}_unknown.fastq.gz",'wt')
output_files[f"{prefix}R2{suffix}_unknown.fastq.gz"] = gzip.open(f"{out_dir}/{prefix}R2{suffix}_unknown.fastq.gz",'wt')

#main loop
while True:
    #as os.listdir returns in arbitrary order, we can't loop over input_list for our filenames
    #instead, reconstruct each filename to take a record from, and loop over the four lines of a FASTQ record
    for i in range(4):
        r1[i] = input_files[f"{prefix}R1{suffix}.fastq.gz"].readline().strip()
        r2[i] = input_files[f"{prefix}R2{suffix}.fastq.gz"].readline().strip()
        r3[i] = input_files[f"{prefix}R3{suffix}.fastq.gz"].readline().strip()
        r4[i] = input_files[f"{prefix}R4{suffix}.fastq.gz"].readline().strip()
    #if we've consumed all FASTQ records, and readline only returns empty strings
    if not r1[0]:
        break
    r3[1] = bioinfo.reverse_complement(r3[1]) # type: ignore
    amend_headers(r1,r2,r3,r4)

    #if indices are unknown, or if they're below a quality cutoff, put into the unknown files
    if (r2[1] not in idx_dict or 
            r3[1] not in idx_dict or 
            bioinfo.qual_score(r2[3]) < 26 or #type: ignore
            bioinfo.qual_score(r3[3]) < 26): # type: ignore
        #increment counter
        bad_idx_ctr += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files[f"{prefix}R1{suffix}_unknown.fastq.gz"].writelines(r1)
        output_files[f"{prefix}R2{suffix}_unknown.fastq.gz"].writelines(r4)
    #else if we have valid but mismatched indices
    elif r2[1] != r3[1]:
        #increment counter
        idx_pair = construct_mismatch(r2[1], r3[1])
        if idx_pair not in mism_ctr:
            mism_ctr[idx_pair] = 1
        else:
            mism_ctr[idx_pair] += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files[f"{prefix}R1{suffix}_mismatched.fastq.gz"].writelines(r1)
        output_files[f"{prefix}R2{suffix}_mismatched.fastq.gz"].writelines(r4)
    #else if we have matching indices
    else: 
        #increment counter
        code: str = idx_dict[r2[1]]
        match_ctr[code] += 1
        #add newlines back in for output
        for i in range(4):
            r1[i] += '\n'
            r4[i] += '\n'
        #write output
        output_files[f"{prefix}R1{suffix}_{code}.fastq.gz"].writelines(r1)
        output_files[f"{prefix}R2{suffix}_{code}.fastq.gz"].writelines(r1)
#close input files
for file in input_files.values():
    file.close()
#close output files
for file in output_files.values():
    file.close()

#TODO: non-file output

for idx_pair in match_ctr.keys():
    if match_ctr[idx_pair] == 0:
        os.remove(f"{out_dir}/{prefix}R1{suffix}_{idx_pair}.fastq.gz")
        os.remove(f"{out_dir}/{prefix}R2{suffix}_{idx_pair}.fastq.gz")

