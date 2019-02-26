# MMG 425 HOMEWORK ASSIGNMENT, SPRING SEMESTER 2019

## LEARNING OBJECTIVES

1. Understand the general process of microbial community analysis using 16S rRNA gene sequence data.  We will use the sequence analysis tools from the Metagenomics Analysis Server (MG-RAST).
2. Familiarize with sequence respositories (where to find and download real data), what "raw" sequence data looks like, and generally how to manipulate these data for analysis using online tools and local worksheets. 
3. Construct charts to visualize and interpret aspects of diversity that are commonly used in community ecology.  
4. Compare bacterial community composition across environmental gradients and hypothesize as to the drivers of the observed patterns.
5.  Reference the literature and resources to build conclusions about the data and suggest next avenues of research.  

## MATERIALS

1. Raw 16S rRNA gene sequence data set from coal mine fire affected soil samples of Centralia.
2. Access to the MG-RAST website. 
3. Spreadsheet and graphics software such as Excel.
4. References about microbes in fire-affected soils and Centralia (available on D2L).

### Getting familiar with the dataset: The Centralia underground coal mine fire!
One of the microbial ecosystems studied by Dr. Shade's lab is Centralia, Pennsylvania. Centralia is the site of an underground coalmine fire that has been burning very slowly along the coal seams since 1962. We study the soils overlying the fire, which are exposed to extreme heat (temperatures have been recorded as hot as 300C! More recently the hottest temperatures are 50-60 C). The objective of our research project is to investigate the impact of an extreme (think: sledgehammer!) and long-term thermal disturbance on soil microbial communities. We use this unusual ecosystem as a model to understand microbial community resilience and stability.  

A neat aspect of the Centralia ecosystem is that as the fire advances along the coal seams, previously impacted soils can recover to ambient temperatures.  This allows for an informative gradient that can be used to figure out how the soil microbiome changes in response to the fire and then responds after the fire subsides locally.  There are also nearby soils that have never been impacted by the fire, but they share very similar geology so that can be used as a "control" for comparison to what is a typical soil microbiome for the soil around Centralia.  

Surface soil cores (overlying the coal seams) were collected in 2014, and classified into three different types according to the fire impact: Fire Affected, Recovered, and Reference. Amplification and high-throughput Illumina sequencing of 16S rRNA gene was conducted from those soil samples to identify the bacterial and archaeal communities. The table below shows the characteristic of the soil and environmental conditions of three sites. The map describes the exact location of the sampling sites in Centralia.

| Sample ID | Soil temperature | Classification | Organic matter | pH | NO3N (ppm) | NH4N (ppm) | Sulfate/Sulfur (ppm) | Arsenic (ppm) | Soil moisture |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| C04 | 13.3 | Recovered | 12.8 | 4.6 | 0.8 | 5 | 23 | 3.6 | 24.2 |
| C13 | 57.4 | Fire Affected | 7.1 | 8 | 4.6 | 1.7 | 28 | 2.58 | 14.7 |
| C17 | 12.1 | Reference | 6.1 | 5.7 | 0.1 | 3.3 | 6 | 1.99 | 13.6 |

![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/centralia_map.png)

The Raw 16S rRNA gene sequences that were obtained from Michigan State's Research Technology Support Facility (RTSF) Genomics Core, here on campus! After data were generated and published, they were made public by submitting them to the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA). 

The Sequence Read Archive (SRA) is an international data repository that stores raw sequencing data and alignment information from high-throughput sequencing platforms (e.g., Illumina MiSeq, HiSeq).  SRA is established under International Nucleotide Sequence Database Collaboration (INSDC). To know more the SRA overview, please see this [link](https://www.ncbi.nlm.nih.gov/sra/docs/). 

Fun fact:  Did you know that tax-payer funded data are often <b>required</b> to be made available to the public and to other researchers?  You can find data from a lot of different microbiome studies on the SRA!

Here, you will also learn how to download the raw sequence data from NCBI, SRA by following the steps below.

### How to fetch your raw sequence data from NCBI Sequence Read Archive (SRA)
1. Go to this [link](https://www.ncbi.nlm.nih.gov/sra/SRP082686).
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/SRA%20Fig%201.png)
There are 54 items and those are raw 16S rRNA gene amplicon sequencing data from Centralia coal mine surface soil consist of 18 different sites and three replicates for each site. 
2. Click the item (for example, No.1 C04D01_GTATGCGCTGTA_L001_R1_001) and you will see the information about the sequence.  
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/SRA%20Fig%202.png)
3. Click the SRR number under "Run" tab (for example: SRR4054183). If you want to download the sequence file, go to "Download" tab and hit "SRA Toolkit" to download SRA Toolkit. The SRA Toolkit is needed to download the sequence and split the paired-ends reads into two fastQ files.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/SRA%20Fig%203.png)
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/SRA%20Fig%204.png)
4. Example of FastQ file downloaded from SRA.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/SRA%20Fig%205.png)

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
**note**: The different format of the original fastQ files and those downloaded from the SRA results in a incompatibility issue for sequence analysis using RDP pipeline. Thus, we will provide you the original fastQ files from Centralia samples (C04, C13, C17) as the input for the next sequence analysis using RDP's Pipeline. There are three paired-end reads (R1 and R2) or six fastQ files compressed into one centralia.tar file that you can download from the D2L.

## GENERAL WORKFLOW FOR SEQUENCE ANALYSIS

<i>Background</i>.
The analysis of 16S rRNA gene sequences is needed to investigate bacterial and archaeal community, their diversity, and structure. The 16S rRNA gene sequence analysis is also commonly conducted to understand the microbial community changes as consequence of the environmental disturbances. 
The figure below will help you to understand the general workflow of microbiome amplicon sequence analysis.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/Workflow_image.png)

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
### 3. Ecological Analyses
1.  There are many programs that can be used to perform diversity analyses and visualization.  A good open-source statistical computing program is R. There are also many custom programs and tools, including QIIME, mothur, RDP-Pipeline, and KBASE.  
2.  Here, we will use MG-RAST and also Excel for basic diversity assessments. 

## TUTORIAL FOR 16S rRNA GENE SEQUENCE ANALYSIS USING MG-RAST PIPELINE
Form into CATME teams and gather together in the class room. Have at least one laptop with wireless connection amongst you, with Excel and a web browser open.
We will provide you three pairs of fastQ files (three tar files) that you can download from D2L.

### Go to this [link](https://www.mg-rast.org/index.html) to register an account on MG-RAST.
Fill the form [here](https://www.mg-rast.org/mgmain.html?mgpage=register).
You will get a notification by email that you have successfully registered an account.

### The MG-RAST Pipeline Steps:
#### Step 1. Initial Processing–Assemble Paired End Reads
1. Download and extract the three tar files on D2L. Each tar file contains two paired fastQ files (forward and reverse). You will use the forward sequence only for each sample ($sequence_name$_L001_R1_001.fastq, e.g C17_06102014_R1_D01_GGCCACGTAGTA_L001_R1_001.fastq).
2. Log in into MG-RAST using your account.
3. Hit "Upload" button as shown below to start upload your forward sequences. You will have three fastQ files for three samples to upload.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig.1.png)
4. Wait until your sequences are succesfully uploaded then hit "Next" button.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig.2.png)
5. Submit your sequences:
1.) Metadata is not required except you want to publish your data. Click "I do not want to supply metadata".
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%203.png)
2.) Put the project name as you want.
3.) Select sequence files.
4.) Choose pipeline options as below.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG-RAST%20Fig%204.png)
5.) Choose the deafault (data stay private) and submit your job.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%205.png)
6. Job Status Monitor will appear and let you know the status of your job.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%206.png)

#### Step 2. Taxonomic Analysis 
1. Go to "analysis" page.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%207.png)
2. Select and add "RDP" to the selected database.
3. Load your project name. Select "name" and your sequences will appear. You also can put your analysis name.
4. Select all of your sequences and load your sequences by hitting the check mark button.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%208.png)
5. Select "Phylum" level in the analysis panel.
![image](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/MG_RAST%20Fig%209.png)
6. Click "Export" button and select "TSV" to download the output and dump it to Excel.
7. Remove the Eukaryota from the Excel file. Examine the output for dominant and rare taxa and answer the question below:
1.  What is the most common taxon at the Phylum level, for each sample? 
2.  What is the rarest taxon at the Phylum level, for each sample?
3.  Make a stacked bar chart. Make a stacked bar for each sample, so you will have three stacked bar in one chart. See the [example](https://github.com/ShadeLab/MMG-425_Assignment/blob/master/Phylum_level_composition.png). Compare:  How does the phylum-level composition vary by soil temperature and fire impact?
4.  Make a hypothesis to explain why particular phyla are found in particular sites.  Hint: use the environmental data provided above, and then use the internet and other references to look up the requirements of the most informative phyla and use that information to inform your hypothesis.  Make sure that you make a statement explaining the phylum distributions for each sample. 

#### Step 3. Making a Species Abundance Distribution to assess evenness.
Use the provided OTU table (there should be 3 columns one for each sample) for this step.  First, sum the total abundance for each OTU.  OTUs are in ROWS, so make a new column of row sums in excel.  Then, sort from high to low.  
1.  From these data, construct a Species Abundance Distribution (also sometimes called a ranked abundance curve). Remember that the x-axis is the ranked taxa from most to least abundant, and the y-axis is the abundance.  Save the chart for the lab report.
2.  Intepret the chart.  What does it tell you about the evenness of the community?

#### Step 4.  Assess the richness for each sample.
Use the provided OTU table for this step.  First, copy the data and paste it into a new sheet. Make a presence-absence table from the data but Finding and Replacing and values > 1 with 1.  Use an IF formula to do so:
=IF(*cell*>1,1,0)
The interpretation of this statement is that if the contents in the cell are greater than 1, replace with 1, and if not, replace with zero.  In the end, the zeros (absenses) will remain zeros.  This transformation removes the relative abundance data and provides presence-absence data, but allows for easy plotting of how many new OTUs are observed with increasing sequencing effort. 
1.  What is the richness for each sample?  Which sample has the highest richness?  Which has the lowest? 
2.  Hypothesize as to what drives the richness pattern that you observe.  What about the environment or ecology would lead to this pattern?

#### Step 5. Make the final report
Make your final report in word document. The final report should be turned in (one per team) on D2L.  The content should include: 
1.  Introduction:  summarize and synthesize your background reading on the Centralia environment and unique expectations of its microbial ecology.  This should be 2-3 paragraphs.  Use the references provided on D2L, but you are welcome to find additional references.  
2.  Results.  Use the headings above for steps 1-4 and answer each question in full sentences. Insert any charts and figures into the report.  
3.  Conclusions and Future Directions.  Based on your exploration of these sequences, what conclusions can you make about the communities that live in Centralia?  What outstanding hypothese do you have?  What additional sampling or experiments could be done to address these hypotheses?  
4.  Annotated references.  Use the numeric alphabetica citation format (the same as for the Microbial Ecosystems Project) and annotated in the same way.  
