# MMG 425 HOMEWORK ASSIGNMENT, SPRING SEMESTER 2019

## LEARNING OBJECTIVE

1. To acquaint MMG 425 participants with the process of microbial community analysis of 16S rRNA gene sequence data using RDP's pipeline and community comparison systems.
2. To familiarize MMG 425 participants with FastQ files and the Ribosomal Database II Project (RDP II) Classifier.
3. To classify communities into bacterial taxa. 
4. To construct: (a) Sampling Curve and (b) Rank Abundance Curve figures from classified 16S rRNA gene sequence files.
5. To compare bacterial community composition of different ecosystems.

## MATERIALS

1. Raw 16S rRNA gene sequence data set from coal mine fire affected soil samples of Centralia.
2. Access to the Ribosomal Database II Project website. 
3. Laptops with wireless connection capability.
4. Spreadsheet and graphics software such as Excel.

### Getting familiar with the dataset: The Centralia underground coal mine fire!

One of the projects in Dr. Shade's lab is Centralia. It is an underground coalmine fire in Pennsylvania that has been burning since 1962. The objective of the project is to investigate the impact of an extreme and long-term thermal disturbance on soil microbial communities. The soil samples were taken and classified into three different types: Fire Affected, Recovered, and Reference. Amplification and sequencing of 16S rRNA gene was conducted from those soil samples to identify the bacterial and archaeal communities. 

Raw 16S rRNA gene sequences that were obtained from the RTSF, MSU sequencing facility, were submitted to the National Center for Biotechnology Information (NCBI) and they are available in NCBI’s Sequence Read Archive (SRA).
The Sequence Read Archive (SRA) is an international bioinformatics database that stores raw sequencing data and alignment information from high-throughput sequencing platforms such as Ilumina MiSeq and it is established under International Nucleotide Sequence Database Collaboration (INSDC). To know more the SRA overview, please see this [link](https://www.ncbi.nlm.nih.gov/sra/docs/)

### How to fetch your raw sequence data from NCBI Sequence Read Archive (SRA)

1. Go to this [link](https://www.ncbi.nlm.nih.gov/sra/SRP082686).
2. There are 54 items and those are raw 16S rRNA gene amplicon sequencing data from Centralia coal mine surface soil consist of 18 different sites and three replicates for each site. Pick three data set, number 1, 8, and 12.
(1) C04_05102014_R1_D01_GTATGCGCTGTA_L001_R1_001 represent Recovered site
(8) C17_06102014_R1_D03_CTAGCGAACATC_L001_R1_001 represent Reference site
(12) C16_06102014_R1_D03_ACGCCACGAATG_L001_R1_001 represent Fire Affected site
3. Click the item and you will see the information about the sequence.
4. Click the accession/SRR number under "Run" tab (for example: SRR4054183). 
5. Go to "Download" tab and hit "SRA Toolkit" to download SRA Toolkit. The SRA Toolkit is needed to download the sequence and split the paired-ends reads into two fastQ files.
6. Download the toolkit version 2.9.2 according to your computer spec.
7. Extract the .tar file in your local computer (for example: sratoolkit.2.9.2-mac64.tar) then put the extracted file in the directory you want.
8. Open the Terminal in your mac (or any command line program that you can install and use for Windows, such as MobaXterm). The app is in the Utilities folder in Applications.
9. Path to the directory where you put the sratoolkit file.
10. Download raw sequence of interest and convert the SRA file into fastQ file by running the two commands below:
```
sratoolkit.2.9.2-mac64/bin/./prefetch <SRRnumber>
sratoolkit.2.9.2-mac64/bin/./fastq-dump --origfmt -I --split-3 <SRRnumber>
```
The first command will give you the SRA file according to the SRR number. The file is located in the ncbi file in your home directory (for example: SRR4054183.sra).
The second command will splits paired reads into two separate fastQ files: *_1.fastq and *_2.fastq.
11. Use these fastQ files as the input for the next sequence analysis using RDP's Pipeline.

### Getting familiar with 16S rRNA sequence data file format:  the FastQ File

<i>Background</i> 
FastQ file is a text file that contains the sequence data generated from the Ilumina sequencing technology. The fastQ file contains four lines as the example of the original FastQ file below. For FastQ file structure, read this [link](https://en.wikipedia.org/wiki/FASTQ_format)
```
@HWI-M03127:41:ACE13:1:1101:15519:1371 1:N:0:GTATGCGCTGTA
TACGAGGGGGGCAAGCGTTGTTCGGAATTATTGGGCGTAAAGGGAGCGTAGGCGGTTCGGTAAGTCACTTGTGAAATCTCTGGGCTCAACTCAGAGTCTGCAAGCGAAACTGCCGGGCTGGAGTATGGGAGAGGTGAGTGGAATTCCTGG
+
BBBBABBBBBBBGFFGGFGGGGHGGEFFHFHHHHHGFGGDHGHGGGGGGGGGHGGGGGGGGHGHHHHHHHHHHGGHGHHHHHHGHHHGHGHHHHHGHHGGGGGGGGGGGGGGGGGGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
FastQ files downloaded from SRA have slightly different header format from the original Ilumina fastQ files because SRA does not store the unique identifiers/read names. The identifiers were replaced with serial numbers as below. You can find the reasons in this [link](https://github.com/ncbi/sra-tools/wiki/Read-Names)
 and [link](https://github.com/ncbi/sra-tools/issues/130)
```
@SRR4054183.1.1 1 length=150
TACGAGGGGGGCAAGCGTTGTTCGGAATTATTGGGCGTAAAGGGAGCGTAGGCGGTTCGGTAAGTCACTTGTGAAATCTCTGGGCTCAACTCAGAGTCTGCAAGCGAAACTGCCGGGCTGGAGTATGGGAGAGGTGAGTGGAATTCCTGG
+SRR4054183.1.1 1 length=150
BBBBABBBBBBBGFFGGFGGGGHGGEFFHFHHHHHGFGGDHGHGGGGGGGGGHGGGGGGGGHGHHHHHHHHHHGGHGHHHHHHGHHHGHGHHHHHGHHGGGGGGGGGGGGGGGGGGGGFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
```
## GENERAL WORKFLOW FOR SEQUENCE ANALYSIS

<i>Background</i> 
The analysis of 16S rRNA gene sequences is needed to investigate bacterial and archaeal community, their diversity, and structure. The 16S rRNA gene sequence analysis is also commonly conducted to understand the microbial community changes as consequence of the environmental disturbances. General workflow for amplicon analysis is divided into four main steps consists of:

1. Data Pretreatment 
Data pretreatment consists of multi-step proccesses including:
(1) <i>Merging paired-end reads</i>. Merging or assembly of paired-end reads generates single fastQ file from forward and reverse fastQ files. A pair is merged by aligning the forward read sequence to the reverse-complement of the reverse read sequence. In the overlap region where both reads cover the same bases, a single letter and Q score is derived from the aligned pair of letters and Q scores for each base. Read this [link](https://drive5.com/usearch/manual9/merge_pair.html) for more information.
(2) <i>Primers and adapters removal</i>. Primers, adapters, or any unwanted type of synthetic sequences should be removed from the sequencing reads, otherwise they will cause sequence contamination.
(3) <i>Quality filtering</i>. Quality filtering is conducted to reduce the number of errors and to avoid adverse effects in the downstream analysis. This can be done by removing reads with poor quality bases (e.g. low quality (Phred) score) and discarding ambiguous/unknown bases (N). 
Quality (Phred) score: the quality score of a base, which is an integer value representing the estimated probability of an error, i.e. that the base is incorrect. See this [link](https://en.wikipedia.org/wiki/Phred_quality_score).
(4) <i>Length trimming</i>. Length trimming is conducted so all of the sequences have the same length and start at the same position to get good OTU sequences. However, length trimming is usually not necessary for overlapping paired-end reads. Trimming for fungal ITS amplicon reads is also not necessary because the nature of ITS gene that have large variations in length.
(5) <i>Denoising & chimera checking</i>. Denoising is removing noise such as sequence errors and chimeras that occurs during PCR amplification or sequencing process.
(6) <i>Eukaryotes contamination removal</i>. Eukaryotes contamination such as chloroplast and mitochondria have to be removed from the amplicon reads.
(7) <i>Rarefaction</i>. Data normalization of quality filtered sequence 

2. OTU table construction 
3. data diversity analysis & visualization (alpha & beta diversity, multivariate analyses, barplot/heatmap etc.), and
4. OTU occupancy & co-occurence analyses.
Data pretreatment and OTU table construction (step 1 & 2) can be conducted using different pipelines/platforms. There are several platforms that can be used for amplicon analysis such as [QIIME](http://qiime.org/tutorials/index.html), [mothur](https://mothur.org/wiki/MiSeq_SOP), [USEARCH](https://www.drive5.com/usearch/manual/uparse_pipeline.html), and the [RDPipeline](http://pyro.cme.msu.edu)(Ribosomal Database Project Pipeline).
Some of those platforms such as QIIME and mothur also can be used to do diversity analyses and some data visualization. However, it is much better and more flexible to do ecological analyses (step 3 & 4) using ecological analysis softwares (vegan, Phyloseq) on R programming language.

















## TUTORIAL FOR 16S rRNA GENE SEQUENCES ANALYSIS USING RDP's PIPELINE

Form into CATME teams and gather together in the class room. Have at least one laptop with wireless connection amongst you, with Excel and a web browser open.

We will use paired-end of 16S rRNA gene sequences that have been downloaded from NCBI, SRA as described previously. There are three pairs of fastQ files (six sequence fastQ files in total) that were compressed into one .tar file. 

### Go to http://rdp.cme.msu.edu/tutorials/workflows/16S_supervised_flow.html to observe the 16S supervised workflow of RDP's pipeline 
To get more info of the RDPipeline, you can go to the [website](http://pyro.cme.msu.edu)

### The RDPipeline Processing Steps:

#### Step 1. Initial Processing (Assemble Paired End Reads)
The initial processing of RDPipeline contains multi-steps process and quality filtering including sorting the raw reads by sample tag, trimming off tag and primer regions, and removing low quality sequences (Cole et al. 2014). It also uses a tool called Assembler to assembly the paired-end reads (Cole et al. 2014).

1.) Use the fastQ files of 16S rRNA gene sequence data that have been downloaded as the input fastQ files
2.) Go to https://pyro.cme.msu.edu/login.spr (you don't need to login or make an account) and hit "TEST DRIVE" button.
3.) Upload the .tar file to the initial processing tool main page: https://pyro.cme.msu.edu/init/form.spr
4.) Fill the initial processing form. The forward primer(s) or reverse primer(s) are not required. 
Job name = use any name
Max number of N's = 0
Min Read Q score = 27
Min sequence length = 250
Max sequence length = 280
5.) Make sure to check the box "Assemble paired end reads". 
6.) Hit "Perform Initial Processing" button and wait for a while.
7.) Check if your job is still running or complete by hitting "my jobs" tab on the right corner of the window.
8.) Download the output file once the job is complete.
9.) Inspect your output data by looking at the output file: http://rdp.cme.msu.edu/tutorials/init_process/RDPtutorial_INITIAL-PROCESS_pe.html

#### Step 2. Taxonomic Classification Using RDP Classifier

1.) Retrieve the FastA files from the D2L course site. 
2.) Open the RDP Classifier at https://pyro.cme.msu.edu/classifier/form.spr and inspect it.
Toggle on the Classifier tab and inspect it. 
3.) “Run” the classifier with each of the three files by attaching the text file as the procedure indicates, and uploading. Wait for the output.  
4.) Examine the output for dominant and rare taxa. 
5.) ANSWER THIS QUESTION: what is the most common taxon at the Phylum level?  
6.) Capture the output in a text file and dump it to Excel.  
8.) Sort the file again by “Phylum” to capture the frequency of matches to Phylum. 
9.) SAVE THE THREE CLASSIFICATIONS IN AN EXCEL FILE USING PHYLUM AS THE IDENTIFIER.  DISCARD ALL OTHER IDENTIFIERS (FAMILY, ORDER, SPECIES, ETC.).

### Step 3:  Making a rank abundance curve to assess evenness
 
1.) Construct a Rank Abundance Curve of Phylum x frequency of match using the histogram graph function in Excel. 
2.) Use only phylum classified results.  
3.) Repeat and compare for the three data sets (C04, C17, C16). 
4.) SAVE THE ABUNDANCE CURVES IN AN EXCEL FILE.

### Step 4:  SAMPLING CURVE: ESTIMATION OF RICHNESS

(1) Construct a Sampling Curve at the genus level using the following procedure. 
(2) Resample each of the FastA files by taking approximately 40 rows of data at a time and classifying each set of 40 rows. 
(3) Use only genus classified results.  Construct a set of samples (for 371 rows, that will be about 9 samples, maybe 10) into a sampling curve by plotting the number of genera that are discovered on the Y axis and the sample number on the X axis.  
(4) ANSWER THESE QUESTIONS:  Does the curve deflect?  Can you estimate Richness from the curve using the asymptotic method? 
(5) SAVE THE FOUR CURVES IN AN EXCEL FILE.

### Make a bar chart of phylum-level composition in excel

### Step 5:  COMMUNITY COMPARISONS using RDP II Libcompare

(1) Using the LIBCOMPARE tool in RDP II, submit pairwise FastQ files for community composition for the four files.  
There will be 6 comparisons.  
The “significance” test is a probability value, thus 6.1EXP-1 is read as P = 0.61 (i.e., not significant).  
A value of 6.1EXP-3 is read as P = 0.0061 (i.e., a significant difference in frequencies of the taxa of interest). 
(2) ANSWER THESE QUESTIONS: Which of the communities are rather similar and which are very different?  Does “organic” differ from “conventional?”  Does “fly gut” differ from “cow feces?”


## MAKE YOUR FINAL REPORT







