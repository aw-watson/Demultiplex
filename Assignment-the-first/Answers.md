# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 | 101 | Phred+33 |
| 1294_S1_L008_R2_001.fastq.gz | index1 | 8 | Phred+33 |
| 1294_S1_L008_R3_001.fastq.gz | index2 | 8 | Phred+33 |
| 1294_S1_L008_R4_001.fastq.gz | read2 | 101 | Phred+33 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**
    
## Part 2
1. Define the problem

-We need to separate FASTQ records based on their indices (stored in R2 and R3, for Illumina output). In this case, we have dual-matched indices, so we want to separate records that have matching indices into separate FASTQ files. 

-After separating these records, we will need to be able to tell which index corresponds to a collection of records. There should be an indicator for this in the name of the output file, as well as a way to double-check for correct separation in the file itself. 

-We will also need to keep records of reads from one direction separate from records of reads from the other direction. 

-For correctness and ease of use, we will also need to separate out records with indices that do not correspond to our known indices or are too low-quality to correctly identify. 

-To measure the extent of index-hopping, we will also need to separate out records with valid but mismatched indices. 

-We will have to track the number of occurrences of all of these types of records, since it would be unfeasible to do so manually. 

2. Describe output

-Up to 48 FASTQ files with names following the format ```1294_S1_L008_R<1|2>_001_<index code>.fastq```.

-2 FASTQ files with names following the format ```1294_S1_L008_R<1|2>_001_unknown.fastq```

-2 FASTQ files with names following the format ```1294_S1_L008_R<1|2>_001_mismatched.fastq```

-All output FASTQ files should append both indices for a record to the header lines of both records in the output.

-1 tab-separated text file named ```mismatched_counts.tsv``` where the first field of each line is the combination of mismatched indices and the second field is a count of how many occurrences of that mismatch were observed in the input.

-1 tab-separated text file named ```matched_counts.tsv``` where the first field of each line is an index code and the second field is a count of how many occurrences of matches of that index pair were observed in the input.

-A printout of form:

```
Number of unknown or low-quality index pairs:    <#>
Total number of mismatched index pairs:    <#>
Total number of matching index pairs:    <#>
```

3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode

**Process:**
1. Create a dictionary of {index sequence:index code} for our known indices.
2. Create a dictionary of {index code: occurrence count} for our known indices, with all values initially 0.
3. Create an empty dictionary to add mismatched pairs to.
4. Initialize a counter for unknown index pairs.
5. Create 4 empty lists of 4 elements to store our records.
6. Create 
7. Open the 4 input files to read from.
8. Open each of the 52 output files to write to.
9. Read one record at a time from each input file (consume 4 lines from each file), storing in the lists.
10. Replace the sequence of the second index with its reverse complement, in place.
11. Add the first index and second index to the headers of both the first and second read.
12. If either index is not in our dictionary of indices OR the average quality score for either index is below a chosen threshold:
      1. Write the record for the first read to the R1 file for unknown reads
      2. Write the record for the second read to the R2 file for unknown reads
      3. Increment a counter for unknown index pairs
      4. Go on to the next record
13. Otherwise, if the indices are not the same:
      1. Write the record for the first read to the R1 file for mismatched reads.
      2. Write the record for the second read to the R2 file for mismatched reads.
      3. Increment a counter (in a dictionary) for mismatched index pairs.
      4. Go on to the next record.
14. Otherwise, we have matching indices:
      1. Write the record for the first read to the R1 file for that index.
      2. Write the record for the second read to the R2 file for that index.
      3. Increment a counter (in a dictionary) for matching index pairs.
      4. Go on to the next record.
15. After processing all records, close all input and output files.
16. Print the number of unknown or low-quality index pairs.
17. Go over the dictionary of mismatched index pairs.
      1. Write each mismatch and the number of occurrences to a tab-separated text file, with one line per unique mismatch.
      2. Print the total number of mismatched index pairs.
18. Go over the dictionary of matching index pairs.
      1. Write each index and the number of occurrences to a tab-separated text file, with one line per index.
      2. If the number of occurrences is 0, delete the R1 and R2 output files corresponding to that index.
      3. Print the total number of matching index pairs.

**End process description**

5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement

Independent Functions:

```python
      def reverse_complement(seq: str) -> str:
          '''Takes a DNA sequence and returns the reverse complement'''
          return rc
      Input: "ATCCCGNNG"
      Expected output: "CNNCGGGAT"
```
      
  ```python
      def amend_headers(R1: list, R2: list, R3: list, R4: list):
          '''Helper function to add index sequences onto header lines of a FASTQ record.
          R2 and R3 are the FASTQ records for the index reads, R1 and R4 are the FASTQ records for the insert reads.'''
          return
      Input: 
             R1: [@, AAA, +, JJJ]
             R2: [@, C, +, J]
             R3: [@, G, +, E]
             R4: [@, TAG, +, JJE]
      Expected end state: 
             R1: [@:C-G, AAA, +, JJJ]
             R2: [@, C, +, J]
             R3: [@, G, +, E]
             R4: [@:C-G, TAG, +, JJE]
  ```

```python
    def qual_score(phred_score: str) -> float:
        """Returns the average Phred score of a string of contiguous Phred+33 scores"""
        return average_phred_score

    Input: "#AAAFJJJ"
    Output: 32.25
```
