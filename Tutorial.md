# MMG 425 HOMEWORK ASSIGNMENT, SPRING SEMESTER 2019

## lEARNING OBJECTIVE
```
1. To acquaint MMG 425 participants with the process of microbial community analysis through manipulation of FastA files with 16S rDNA nucleotide sequence data and the RDP II classifier and community comparison systems.
2. To familiarize MMG 425 participants with FastA files and the Ribosomal Database II Project (RDP II) Classifier.
3. To classify communities to bacterial taxa. 
4. To construct:
(a) Sampling Curve and 
(b) Rank Abundance Curve 
figures from classified FastA 16S rDNA nucleotide sequence files.
5. To compare bacterial community composition of different ecosystems.
```
## MATERIALS
```
1. Four FastA files of 16S rDNA nucleotide sequence data from a comparative study done on dairy farms in Wisconsin: (a) feces from dairy cows on organic farms; (b) feces from cows on conventional farms; (c) guts of flies collected on organic dairy farms; (d) guts of flies collected on conventional dairy farms.  
2. Access to the Ribosomal Database II Project website. 
3. Laptops with wireless connection capability.
4. Spreadsheet and graphics software such as Excel.
```
### FastA Files
```
For FastA file structure, read https://submit.ncbi.nlm.nih.gov/genbank/help/ and http://en.wikipedia.org/wiki/FASTA_format

FastA files are simple text files containing a descriptor of a nucleotide or amino acid sequence set apart by a > sign and a name of the sequence, followed on the next line by the sequence itself.  There must be carriage returns at the end of each line, and each line may not be more than 80 characters in length.

Example of three 16S rDNA nucleotide sequences in the FastA format:

>FC102.01
GGGGACCTTCGGGCCTTGCGCTATCAGATGAGCCTAGGTCGGATTAGCTAGTTGGTGAGGTAATGGCTCCCCAAGGCTACGATCCGTAACTGGTCTGAGAGGATGATCAGTCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCAGGCCGCGTGTGTGAAGAAGGTCTTCGGATTGTAAAGCAAGGAAGTTGGGAGGAAGGGC
>FC102.02
CACGGATAACATACCGAAAGGTATGCTAATACGAGATAATATGCTTTTATCGCACGGTAGAAGTATCAAAGCTCCGGCGGTACAGGATGGACCCGCGTCTGATTAGCTAGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCAGTAGCCAACCTGAGAGGGTGATCGGCCACATTGGAACTGAGACACGGTCCAAACTCTTACGGGAGGCAGCAGTGGGGAATATTTTAAAAGGGGCGAAAGCCTGATGCAGCA
>FC102.03
CCCTAGAGTGGGGGATAACGTAGCGAAAGTTACGCTAATACCGCATACGATCTAAGGGTGAAAGTGGGGGATCGCAAGACCTCATGCTCGTGGAGCGGCCAATATCTGATTAGCTAGTTGGTAGGGTAAAAGCCTACCAAGGCATCGATCAGTAGCTGGTCTGAGAGGACGACCAGCCACACTGGAGCTGAGACACGGTCCACACTCCTACGGGAGGCAGCAGTGGGGAATTTTGGACAATGGGCGAAAGCCTGA
```
## TUTORIAL
```
(1) Form into CATME teams and gather together in the class room.  Have at least one laptop with wireless connection amongst you, with Excel and a web browser open.

(2) CLASSIFICATION.  Retrieve the FastA files from the D2L course site.  Open the RDP II website at http://rdp.cme.msu.edu/ and inspect it.  Toggle on the Classifier tab and inspect it.  “Run” the classifier with each of the four files by attaching the text file as the procedure indicates, and uploading.  Wait for the output.  Examine the output for dominant and rare taxa.  ANSWER THIS QUESTION: what is the most common taxon at the genus level?  Capture the output in a text file and dump it to Excel.  Sort the file to remove the “zero” values.  Sort the file again by “genus” to capture the frequency of matches to genus.  SAVE THE FOUR CLASSIFICATIONS IN AN EXCEL FILE USING GENUS AS THE IDENTIFIER.  DISCARD ALL OTHER IDENTIFIERS (FAMILY, ORDER, SPECIES, ETC.).

(3) RANK ABUNDANCE CURVE: EVENNESS.  Construct a Rank Abundance Curve of genus x frequency of match using the histogram graph function in Excel.  Use only genus classified results.  Repeat and compare for the remaining three FastA files.  SAVE THE FOUR CURVES IN AN EXCEL FILE.

(4) SAMPLING CURVE: ESTIMATION OF RICHNESS.  Construct a Sampling Curve at the genus level using the following procedure.  Resample each of the FastA files by taking approximately 40 rows of data at a time and classifying each set of 40 rows.  Use only genus classified results.  Construct a set of samples (for 371 rows, that will be about 9 samples, maybe 10) into a sampling curve by plotting the number of genera that are discovered on the Y axis and the sample number on the X axis.  ANSWER THESE QUESTIONS:  Does the curve deflect?  Can you estimate Richness from the curve using the asymptotic method?  SAVE THE FOUR CURVES IN AN EXCEL FILE.

(5) COMMUNITY COMPARISONS.  Using the LIBCOMPARE tool in RDP II, submit pairwise FastA files for community composition for the four files.  There will be 6 comparisons.  The “significance” test is a probability value, thus 6.1EXP-1 is read as P = 0.61 (i.e., not significant).  A value of 6.1EXP-3 is read as P = 0.0061 (i.e., a significant difference in frequencies of the taxa of interest).  ANSWER THESE QUESTIONS: Which of the communities are rather similar and which are very different?  Does “organic” differ from “conventional?”  Does “fly gut” differ from “cow feces?”
```









