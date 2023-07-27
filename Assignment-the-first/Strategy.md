Description of de-multiplexing process:

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
          R2 is the record for the index read corresponding to R1's sequence read.
          R4 is the record for the index read corresponding to R3's sequence read.'''
          
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

Variables:

1. A list of input filenames, named ```in_list```

   ```
   1294_S1_L008_R1_001.fastq.gz
   1294_S1_L008_R2_001.fastq.gz
   1294_S1_L008_R3_001.fastq.gz
   1294_S1_L008_R4_001.fastq.gz
   ```

2. A list of output filenames, named ```out_list```
   ```
   1294_S1_L008_R1_001_nonmatching.fastq
   1294_S1_L008_R2_001_nonmatching.fastq
   1294_S1_L008_R1_001_unknown.fastq
   1294_S1_L008_R2_001_unknown.fastq
   1294_S1_L008_R1_001_B1.fastq
   1294_S1_L008_R2_001_B1.fastq
   1294_S1_L008_R1_001_B9.fastq
   1294_S1_L008_R2_001_B9.fastq
   .
   .
   .
   1294_S1_L008_R<1|2>_001_<Index ID>.fastq
   ```

3. A dictionary mapping index sequences to their respective identifiers, named ```idx_codes```

4. A dictionary mapping all possible matching index codes to observed counts, named ```mtch_ctr```

5. A dictionary mapping all possible mismatched pairs of index codes to observed counts, named ```msmtch_ctr```. Pairs of mismatched codes should be in alphabetical order to avoid duplicates.

6. An accumulator (int) counting the number of read pairs with unknown indexes, named ```unk_ctr```.

**Process:**

```
for each name in in_list:
  -open the file in read mode
for each name in out_list:
  -open the file in write mode

while the input files have lines:
  -set unk_ctr, all values in mtch_ctr, and all values in msmtch_ctr to 0
  -consume 4 lines from each input file, and store them in lists: R1, R2, R3, R4
  -replace the second index read (in R3) with its reverse complement, storing it back into the list
  -pass lists to amend_headers

  if the first index read is not in idx_codes
    OR the second index read is not in idx_codes
    OR the first index read's average quality score is below a cutoff
    OR the second index read's average quality score is below a cutoff:

    -write R1 to 1294_S1_L008_R1_001_unknown.fastq
    -write R4 to 1294_S1_L008_R2_001_unknown.fastq
    -increment unk_ctr
    -continue to next iteration

  if the two index reads are not the same sequence:
    -write R1 to 1294_S1_L008_R1_001_nonmatching.fastq
    -write R4 to 1294_S1_L008_R2_001_nonmatching.fastq
    -increment the appropriate entry in msmtch_ctr
    -continue to next iteration

  #if execution reaches this line, we have matching indexes
  increment the appropriate entry in mtch_ctr
  write R1 to 1294_S1_L008_R1_001_<index code>.fastq
  write R4 to 1294_S1_L008_R2_001_<index code>.fastq

for each name in in_list:
  -close file
for each name in out_list:
  -close file

for each entry in mtch_ctr:
  if the observed count is 0:
    -delete 1294_S1_L008_R1_001_<index code>.fastq
    -delete 1294_S1_L008_R2_001_<index code>.fastq

-print number of index pairs with unknown indices
for each entry in mtch_ctr:
  -print the number of matching index pairs corresponding to that code
for each entry in msmtch_ctr:
  -print the number of index pairs with that mismatch #if possible, format this as a triangular table
```
