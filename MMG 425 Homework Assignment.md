# MMG 425 HOMEWORK ASSIGNMENT, SPRING SEMESTER 2019

## LEARNING OBJECTIVE

1. Understand the general process of microbial community analysis using 16S rRNA gene sequence data.  We will use the sequence analysis tools from the Ribosomal Database Project (RDP).
2. Familiarize with sequence respositories (where to find and download real data), what "raw" sequence data looks like, and generally how to manipulate these data for analysis using online tools and local worksheets. 
3. Construct charts to visualize and interpret aspects of diversity that are commonly used in community ecology.  
4. Compare bacterial community composition across environmental gradients and hypothesize as to the drivers of the observed patterns.
5.  Reference the literature and resources to build conclusions about the data and suggest next avenues of research.  

## MATERIALS

1. Raw 16S rRNA gene sequence data set from coal mine fire affected soil samples of Centralia.
2. Access to the Ribosomal Database II Project website. 
3. Laptops with wireless connection capability.
4. Spreadsheet and graphics software such as Excel.

### Getting familiar with the dataset: The Centralia underground coal mine fire!
One of the projects in Dr. Shade's lab is Centralia. It is an underground coalmine fire in Pennsylvania that has been burning since 1962. The objective of the project is to investigate the impact of an extreme and long-term thermal disturbance on soil microbial communities. The soil samples were taken and classified into three different types: Fire Affected, Recovered, and Reference. Amplification and sequencing of 16S rRNA gene was conducted from those soil samples to identify the bacterial and archaeal communities. The table below shows the characteristic of the soil and environmental condition of the three sites where the three samples were taken. The map describes the exact location of the sampling sites in Centralia.

| Sample ID | Soil temperature | Classification | Organic matter | pH | NO3N (ppm) | NH4N (ppm) | Sulfate/Sulfur (ppm) | Arsenic (ppm) | Soil moisture |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C04 | 13.3 | Recovered | 12.8 | 4.6 | 0.8 | 5 | 23 | 3.6 | 24.2 |
| C13 | 57.4 | Fire Affected | 7.1 | 8 | 4.6 | 1.7 | 28 | 2.58 | 14.7 |
| C17 | 12.1 | Reference | 6.1 | 5.7 | 0.1 | 3.3 | 6 | 1.99 | 13.6 |

![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/centralia_map.png)

The Raw 16S rRNA gene sequences that were obtained from the Research Technology Support Activity (RTSF), MSU sequencing facility, were submitted to the National Center for Biotechnology Information (NCBI) and they are available in NCBI’s Sequence Read Archive (SRA).
The Sequence Read Archive (SRA) is an international bioinformatics database that stores raw sequencing data and alignment information from high-throughput sequencing platforms such as Ilumina MiSeq and it is established under International Nucleotide Sequence Database Collaboration (INSDC). To know more the SRA overview, please see this [link](https://www.ncbi.nlm.nih.gov/sra/docs/). 
Here, you will also learn how to download the raw sequence data from NCBI, SRA by following the steps below.
### How to fetch your raw sequence data from NCBI Sequence Read Archive (SRA)
1. Go to this [link](https://www.ncbi.nlm.nih.gov/sra/SRP082686).
2. There are 54 items and those are raw 16S rRNA gene amplicon sequencing data from Centralia coal mine surface soil consist of 18 different sites and three replicates for each site. Pick the data that you want to download. For example, number 1. [C04D01_GTATGCGCTGTA_L001_R1_001](https://www.ncbi.nlm.nih.gov/sra/SRX2043754[accn]).
3. Click the item and you will see the information about the sequence.
4. Click the SRR number under "Run" tab (for example: SRR4054183). 
5. Go to "Download" tab and hit "SRA Toolkit" to download SRA Toolkit. The SRA Toolkit is needed to download the sequence and split the paired-ends reads into two fastQ files.
6. Download the toolkit version 2.9.2 according to your computer spec.
7. Extract the .tar file in your local computer (for example: sratoolkit.2.9.2-mac64.tar) then put the extracted file in the directory you want.
8. Open the Terminal in your mac (or any command line program that you can install and use for Windows, such as MobaXterm). The app is in the Utilities folder in Applications.
9. Path to the directory where you put the sratoolkit file.
10. Download raw sequence of interest and convert the SRA file into fastQ file by running the two commands below:
```
sratoolkit.2.9.2-mac64/bin/./prefetch <SRRnumber>
sratoolkit.2.9.2-mac64/bin/./fastq-dump --skip technical -I --split-3 <SRRnumber>
```
The first command will give you the SRA file according to the SRR number. The file is located in the ncbi file in your home directory (for example: SRR4054183.sra).
The second command will splits paired reads into two separate fastQ files: *_1.fastq and *_2.fastq.

### Getting familiar with 16S rRNA sequence data file format: the FastQ File
<i>Background</i>. 
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
**note**: The difference format of the original fastQ files and those downloaded from the SRA results in the incompatibility issue for sequence analysis using RDP pipeline. Thus, we will provide you the original fastQ files from Centralia samples (C04, C13, C17) as the input for the next sequence analysis using RDP's Pipeline. There are three paired-end reads (R1 and R2) or six fastQ files compressed into one .tar file that you can download from the D2L.

## GENERAL WORKFLOW FOR SEQUENCE ANALYSIS

<i>Background</i>.
The analysis of 16S rRNA gene sequences is needed to investigate bacterial and archaeal community, their diversity, and structure. The 16S rRNA gene sequence analysis is also commonly conducted to understand the microbial community changes as consequence of the environmental disturbances. General workflow for amplicon analysis is divided into three main steps consists of:
### 1. Data Pretreatment 
Data pretreatment consists of multi-step proccesses including:
1. <i>Merging paired-end reads</i>. 
Merging or assembly of paired-end reads generates single fastQ file from forward and reverse fastQ files. A pair is merged by aligning the forward read sequence to the reverse-complement of the reverse read sequence. In the overlap region where both reads cover the same bases, a single letter and Q score is derived from the aligned pair of letters and Q scores for each base. Read this [link](https://drive5.com/usearch/manual9/merge_pair.html) for more information.
2. <i>Primers and adapters removal</i>.
Primers, adapters, or any unwanted type of synthetic sequences should be removed from the sequencing reads, otherwise they will cause sequence contamination.
3. <i>Quality filtering</i>. 
Quality filtering is conducted to reduce the number of errors and to avoid adverse effects in the downstream analysis. This can be done by removing reads with poor quality bases (e.g. low quality (Phred) score) and discarding ambiguous/unknown bases (N). Quality (Phred) score: the quality score of a base, which is an integer value representing the estimated probability of an error, i.e. that the base is incorrect. See this [link](https://en.wikipedia.org/wiki/Phred_quality_score).
4. <i>Length trimming</i>. 
Length trimming is conducted so all of the sequences have the same length and start at the same position to get good OTU sequences. However, length trimming is usually not necessary for overlapping paired-end reads. Trimming for fungal ITS amplicon reads is also not necessary because the nature of ITS gene that have large variations in length.
5. <i>Denoising & chimera checking</i>. 
Noise such as sequence errors and PCR chimeras that occurs during PCR amplification or sequencing process should be eliminated from the amplicon reads because they can cause inflation of microbial diversity estimation.
6. <i>Eukaryotes contamination removal</i>.
Eukaryotes contamination such as chloroplast and mitochondria have to be removed from the amplicon reads.
7. <i>Rarefaction</i>.
Data normalization is needed because the sequencing depth (the number of clean or quality filtered sequences) across samples can significantly differ. Uneven sequencing depth can affect microbial diversity estimation in one sample (alpha diversity). It also can affect microbial diversity comparison among samples.
### 2. OTU table construction 
1. <i>OTU picking</i>.
Operational taxonomic unit (OTU) is cluster or group of similar sequence variants that represent a taxonomic unit of bacteria/archaea species or genus. Read more about OTU in this [link](https://en.wikipedia.org/wiki/Operational_taxonomic_unit). The quality filtered sequences are clustered into OTUs based on the sequence identity cutoffs of 16S rRNA gene of 97, 98, or 99 %. There are three different strategies of OTU clustering namely de novo, closed-reference, and open-reference. Read this [link](http://qiime.org/tutorials/otu_picking.html) for more information.
2. <i>Taxonomy assignment (classifying)</i>.
Taxonomic classification of OTU representatives into seven levels (Domain, Phylum, Class, Order, Family, Genus, Species) can be conducted by alignment against a 16S rRNA gene reference database such as RDP, SILVA database, GreenGenes.
3.  We will provide an OTU table for you for this exercise.
### 3. Ecological Analyses (alpha & beta diversity, multivariate analysis, data visualization, etc.)
Data pretreatment and OTU table construction can be conducted using different pipelines/platforms. There are several platforms that can be used for amplicon analysis such as [QIIME](http://qiime.org/tutorials/index.html), [mothur](https://mothur.org/wiki/MiSeq_SOP), [USEARCH](https://www.drive5.com/usearch/manual/uparse_pipeline.html), and the [RDPipeline](http://pyro.cme.msu.edu)(Ribosomal Database Project Pipeline). Some of those platforms such as QIIME and mothur also can be used to do diversity analyses and some data visualization. However, it is much better and more flexible to do ecological analyses (e.g. microbial diversity, occupancy) using ecological packages (vegan, Phyloseq, microbiome) on R. The figure below will help you to understand the workflow of sequence analysis in general.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/Workflow_image.png)

## TUTORIAL FOR 16S rRNA GENE SEQUENCES ANALYSIS USING RDP's PIPELINE
Form into CATME teams and gather together in the class room. Have at least one laptop with wireless connection amongst you, with Excel and a web browser open.
We will use paired-end of 16S rRNA gene sequences that have been downloaded from NCBI, SRA as described previously. There are three pairs of fastQ files (six sequence fastQ files in total) that were compressed into one .tar file. 

### Go to this [link](http://rdp.cme.msu.edu/tutorials/workflows/16S_supervised_flow.html) to observe the 16S supervised workflow of RDP's pipeline 
To get more info of the RDPipeline, you can go to the [website](http://pyro.cme.msu.edu)

### The RDPipeline Processing Steps:
#### Step 1. Initial Processing–Assemble Paired End Reads
The initial processing of RDPipeline contains multi-steps process and quality filtering including sorting the raw reads by sample tag, trimming off tag and primer regions, and removing low quality sequences (Cole et al. 2014). It also uses a tool called Assembler to assembly the paired-end reads (Cole et al. 2014). Go to this [link]((http://rdp.cme.msu.edu/tutorials/init_process/RDPtutorial_INITIAL-PROCESS_pe.html) to get you familiar with the tool.
1. Use the fastQ files of 16S rRNA gene sequence data that have been downloaded from D2L as the input fastQ files.
2. Go to this [link](https://pyro.cme.msu.edu/init/form.spr).
3. Upload the .tar file to the initial processing tool [main page](https://pyro.cme.msu.edu/init/form.spr).
4. Fill the initial processing form. The taq file and forward primer(s) or reverse primer(s) are not required. 
Job name = use any name that you want.
Max number of N's = 0 (number of ambiguous base allowed).
Min Read Q score = 27 (minimum Phred score used).
Min sequence length = 250.
Max sequence length = 280.
5. Make sure to check the box "Assemble paired end reads". 
6. Hit "Perform Initial Processing" button and wait for a while (be patient..).
7. Check if your job is still running or complete by hitting "my jobs" tab on the right corner of the window.
8. Download the output file once the job is complete.
9. Inspect your output data by looking at the output file [example](http://rdp.cme.msu.edu/tutorials/init_process/RDPtutorial_INITIAL-PROCESS_pe.html).
10. The output file will contain two folders/directories. First folder "assembled_paired_end_sequences" contains the assembled paired sequences in FASTQ format and the assembled analysis results in (FASTQ.LOG) directly output from the Assembler (first stage), one for each pair of input FASTAQ files. Second, the results of stage 2 that are combined in folder called "NoTag". The *NoTag_trimmed.fastq* file is trimmed sequences that are ready for downstream analysis, such as classification using the RDP Classifier.

#### Step 2. Taxonomic Classification Using RDP Classifier
1. Open the [RDP Classifier](http://rdp.cme.msu.edu/tutorials/classifier/classifer_cover_page.html) and inspect it. We will not run RDP Classifier interactively because our data are more than 50 MB. Instead, we will run RDP Classifier on RDPipeline (number 2).
2. Open this [link](http://rdp.cme.msu.edu/tutorials/classifier/RDPtutorial_CLASSIFIER.html) to get you familiar with the tool, input, and output files. 
3. Upload the *NoTag_trimmed.fastq* file in this [link](https://pyro.cme.msu.edu/classifier/form.spr), select "fixrank" and wait for the output. 
4. Examine the output for dominant and rare taxa. Capture the output in a text file and dump it to Excel. Sort the file again by “Phylum” to capture the frequency of matches to Phylum. SAVE THE THREE CLASSIFICATIONS IN AN EXCEL FILE USING PHYLUM AS THE IDENTIFIER AND ANSWER THESE QUESTIONS BELOW:
1.  What is the most common taxon at the Phylum level, for each sample? 
2.  What is the rarest taxon at the Phylum level, for each sample?
3.  Make a stacked bar chart. Make a stacked bar for each sample, so you will have three stacked bar in one chart. See the [example](). Compare:  How does the phylum-level composition vary by soil temperature and fire impact?
4.  Make a hypothesis to explain why particular phyla are found in particular sites.  Hint: use the environmental data provided above, and then use the internet and other references to look up the requirements of the most abundant phyla and use that information to inform your hypothesis.  Make sure that you make a statement explaining the phylum distributions for each sample. 

### Step 3. Making a Species Abundance Distribution to assess evenness.
Use the provided OTU table (there should be 3 columns one for each sample) for this step.  First, sum the total abundance for each OTU.  OTUs are in ROWS, so make a new column of row sums in excel.  Then, sort from high to low.  
1.  From these data, construct a Species Abundance Distribution (also sometimes called a ranked abundance curve). Remember that the x-axis is the ranked taxa from most to least abundant, and the y-axis is the abundance.  Save the chart for the lab report.
2.  Intepret the chart.  What does it tell you about the evenness of the community?

### Step 4.  Assess the richness for each sample.
Use the provided OTU table for this step.  First, copy the data and paste it into a new sheet. Make a presence-absence table from the data but Finding and Replacing and values > 1 with 1.  Us an IF formula to do so:
=IF(*cell*>1,1,0)
The interpretation of this statement is that if the contents in the cell are greater than 1, replace with 1, and if not, replace with zero.  In the end, the zeros (absenses) will remain zeros.  This transformation removes the relative abundance data and provides presence-absence data, but allows for easy plotting of how many new OTUs are observed with increasing sequencing effort. 
1.  What is the richness for each sample?  Which sample has the highest richness?  Which has the lowest? 
2.  Hypothesize as to what drives the richness pattern that you observe.  What about the environment or ecology would lead to this pattern?

### Step 5. Make the final report
Make your final report in word document. The final report should be turned in (one per team) on D2L.  The content should include: 
1.  Introduction:  summarize and synthesize your background reading on the Centralia environment and unique expectations of its microbial ecology.  THis should be 2-3 paragraphs.  Use the references provided on D2L, but you are welcome to find additional references.  
2.  Results.  Use the headings above for steps 1-4 and answer each question in full sentences. Insert any charts and figures into the report.  
3.  Conclusions and Future Directions.  Based on your exploration of these sequences, what conclusions can you make about the communities that live in Centralia?  What outstanding hypothese do you have?  What additional sampling or experiments could be done to address these hypotheses?  
4.  Annotated references.  Use the numeric alphabetica citation format (the same as for the Microbial Ecosystems Project) and annotated in the same way.  







