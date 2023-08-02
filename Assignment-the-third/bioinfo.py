#!/usr/bin/env python
# Author: André Watson, modified from template.


'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.'''

__version__ = "0.5"

DNA_bases = set('ACTGNactgn')
RNA_bases = set('AUCGNaucgn')
DNA_rc_map = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
RNA_rc_map = {'A':'U', 'U':'A', 'G':'C', 'C':'G', 'N':'N'}

def reverse_complement(sequence: str, RNAflag: bool = False) -> str:
    '''Returns the reverse complement of a DNA or RNA sequence'''
    assert validate_base_seq(sequence, RNAflag)
    if RNAflag:
        mapping = str.maketrans(RNA_rc_map)
        sequence = sequence.upper().translate(mapping)
        return sequence[::-1]
    else:
        mapping = str.maketrans(DNA_rc_map)
        sequence = sequence.upper().translate(mapping)
        return sequence[::-1]

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Returns the average Phred score of a string of contiguous Phred scores"""
    sum: int = 0
    for qcall in phred_score:
        sum += convert_phred(qcall)
    return sum/len(phred_score)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(DNA: str):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters - are you sure you used a DNA sequence?"
    
    DNA = DNA.upper()
    return (DNA.count("G")+DNA.count("C"))/len(DNA)

def oneline_fasta(infile: str, outfile: str):
    """Takes a FASTA file (infile) and puts all sequence lines on one line. Writes the result to outfile"""
    with open(infile, 'rt') as rf, open(outfile, 'wt') as wf:
        seq: str = ''
        first_record: bool = True
        while True:
            line = rf.readline().strip()
            if not line:
                break
            if line.startswith(">"):
                if first_record:
                    wf.write(line + "\n")
                    first_record = False
                else:
                    wf.write(seq + "\n")
                    seq = ""
                    wf.write(line + "\n")
            else:
                seq += line
        wf.write(seq)

def calc_median(in_list : list) -> float:
    '''Returns the median of a sorted list of numbers. Undefined on unsorted lists.'''
    if len(in_list) % 2 == 0:
        return (in_list[len(in_list)//2] + in_list[(len(in_list)//2) - 1])/2
    else:
        return in_list[len(in_list)//2]



if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")


    assert validate_base_seq("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", False) #basic operation
    assert validate_base_seq("AUuuuca", True) #mixed case, RNA with uracil
    assert not validate_base_seq("FGGGFGGGGGGGGCF", False) #failure (DNA)
    assert not validate_base_seq("UUUUUUUUnNNNuL", True) #failure (RNA)
    assert validate_base_seq("GGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGCGGGG", True) #RNA without U
    assert validate_base_seq("GGGGGGGGGGGGGGGGGGGNGGGGnGGGGGGGGGG", True) #RNA with N
    assert validate_base_seq("AAAAAAAANNNNNNNNNNNNNnnnAAAA", False) #DNA with N
    print("validate_base_seq appears to be working properly.")
    
    assert gc_content('GGGGGGGGGGGGGGCCCCC') == 1
    assert gc_content('AAAAAAAAAATTTT') == 0
    assert gc_content('GCATGCATCGAT') == 0.5
    assert gc_content('ggccaaatt') == 4.0/9.0
    print("gc_content appears to be working properly.")

    assert calc_median([0,1,2]) == 1 #basic operation
    assert calc_median([7,7,7,7,7]) == 7 #odd list with same entries
    assert calc_median([3,3,6,9,9,9]) == 15/2 #odd list with different entries
    assert calc_median([4.0, 30000.34, 12.1]) == 30000.34 #floats
    assert calc_median([-3, 80.5, 99, 200]) == 179.5/2 #mixed float/int, negative numbers
    print("calc_median appears to be working properly.")

    assert reverse_complement("ACTGGT") == "ACCAGT"
    assert reverse_complement("NNNTCGA") == "TCGANNN"
    assert reverse_complement("ACTGGT", False) == "ACCAGT"
    assert reverse_complement("NNNTCGA", False) == "TCGANNN"
    assert reverse_complement("ACUGGU", True) == "ACCAGU"
    assert reverse_complement("NNNUCGA", True) == "UCGANNN"
    print("reverse_complement appears to be working properly.")