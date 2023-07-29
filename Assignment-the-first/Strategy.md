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

**Process:**
1. Open the 4 input files to read from.
2. Open each of the 52 output files to write to.
3. Read one record at a time from each input file, storing in appropriate lists.
4. Replace the sequence of the second index with its reverse complement.
5. Add the first index and second index to the headers of both the first and second read.
6. If either index is not in our list of indices OR the average quality score for either index is below a chosen threshold:
      1. Write the record for the first read to 
      2. Write the record for the second read to
      3. Increment a counter for unknown index pairs
      4. Go on to the next record
7. Otherwise, if the indices are not the same:
      1. Write the record for the first read to 
      2. Write the record for the second read to
      3. Increment a counter (in a dictionary) for mismatched index pairs
      4. Go on to the next record
8. Otherwise, we have matching indices:
      1. Write the record for the first read to
      2. Write the record for the second read to
      3. Increment a counter (in a dictionary) for matching index pairs
      4. Go on to the next record
9. After processing all records, close all input and output files.
10. Output
