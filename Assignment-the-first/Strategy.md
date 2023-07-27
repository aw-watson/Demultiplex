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
1. Open the 4 input files to read from.
2. Open each of the 52 output files to write to.
3. Read one record at a time from each input file, storing in appropriate lists.
4. Replace the sequence of the second index with its reverse complement.
5. Add the first index and second index to the headers of both the first and second read.
6. If either index is not in our list of indices OR the average quality score for either index is below a chosen threshold:
      1. Write the record for the first read to 
      2. Write the record for the second read to
      3. Increment
      4. Go on to the next record
7. Otherwise, if the indices are not the same:
      1. Write the record for the first read to 
      2. Write the record for the second read to
      3. Increment
      4. Go on to the next record
8. Otherwise, we have matching indices:
      1. Write the record for the first read to
      2. Write the record for the second read to
      3. Increment
      4. Go on to the next record
9. After processing all records, close all input and output files.
10. Output
