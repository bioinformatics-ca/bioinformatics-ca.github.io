---
layout: post2
permalink: /informatics_for_high-throughput_data_sequencing_2016_module7_lab/
title: Informatics for High-Throughput Sequencing Data 2016 Module 7 lab
header1: Informatics for High-Throughput Sequencing Data 2016
header2: Module 7 lab
image: CBW_High-throughput_icon.jpg
---
# CBW HT-seq Module 7 - Galaxy lab
This lab was created by Florence Cavalli

## Table of contents

 [Introduction](#introduction)   
 1. [Load the data](#load)   
 2. [Create single interval](#interval)  
 3. [Check the quality](#quality) 
 4. [Convert the FASTQ quality format](#convert)   
 5. [Trim the read and remove adapter sequence with Trimmomoatic](#trim)   
 6. [Align the reads with BWA-mem](#align)   
 7. [Sort the sam/bam](#sort)   
 8. [Convert bam to sam](#convertBam)   
 9. [Indel realignment](#indelRealign)   
 10. [FixMates](#fixmates)   
 11. [Mark duplicates](#markdup)   
 12. [Base recalibration](#recalibration)   
 13. [Extract Metrics](#extracmetrics)   
 14. [Call SNPs](#callsnp)   
 15. [Filter the variants](#filtervariant)   
 16. [Annotate the variants](#annotate)   
 [The "History options" menu](#menu)   
 [Workflow](#workflow)


## Introduction
<a name="introduction"></a>

In this practical we will use Galaxy to perform the analysis you ran in Module 2 and Module 5.
We are therefore using the dataset from patient NA12878 (reads extracted from region chr1:17704860-18004860)

We will use the galaxy website: https://usegalaxy.org/   
Start by loging on the website with your username

**The paramaters that need to be set or changed are specified in the text below (and in the images), keep the rest of the parameters with their default values**

** indicates the tool to use for each step. You need to select it from the left column   

![tool](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_tools.png) 

#### 1) Load the data 
<a name="data"></a>
** Get Data/Upload File

“Paste/Fetch Data” box:   

- NA12878_CBW_chr1_R1.fastq.gz   
http://cbw##.dyndns.info/HTSeq_module2/raw_reads/NA12878/NA12878_CBW_chr1_R1.fastq.gz

- NA12878_CBW_chr1_R2.fastq.gz   
http://cbw##.dyndns.info/HTSeq_module2/raw_reads/NA12878/NA12878_CBW_chr1_R2.fastq.gz

- hg19_chr1.fa   
http://cbw##.dyndns.info/Module7/hg19_chr1.fa

- dbSNP_135_chr1.vcf.gz   
http://cbw##.dyndns.info/HTSeq_module2/reference/dbSNP_135_chr1.vcf.gz

replace ## by your student id

![data](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_Load_data.png) 

When the data is loaded and/or a task has finished to run, the file created is green. Grey indicates that the task is queueing, yellow that the task is beeing processed and red that an error occurred.   
For example after loading the data, you should see the following in the history column on the right

File check:

![file1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_1.png) 

**NOTE**: Numbers and paths of the files may vary with usage

**NOTE**: The icons on the right of each file allow you to "Poke an eye", edit the attribute or delete the file


#### 2) Create single interval  
<a name="interval"></a>

We create an interval corresponding to the region of interest (we extracted reads from this region for this example dataset) to be used later on by different tools

** Text manipulation/Create Single Interval
 
 - Chromosome: chr1   
 - Start position: 17704860  
 - End position: 18004860    
 - Name: chr1 interval

![createInterval](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_create_interval.png) 

File check:

![file2bis](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_2bis.png) 


#### All the following steps are detailed in the Module 2 practical

#### 3) Check the quality
<a name="quality"></a>

BVAtools is not available on the galaxy website   
You can use ** NGS: QC and manipulation/FastQC Read Quality reports instead   
Since the tool in not avaible we will not run it but this is a **very important step of the analysis**, don't skip it!     
You can use FastQC at several points in your analysis on fastq, bam or sam files


#### 4) Convert the FASTQ quality format
<a name="convert"></a>

** NGS: QC and manipulation/FASTAQ groomer convert between various FASTQ quality formats   
 - File to groom : NA12878_CBW_chr1_R1.fastq

![groomer](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_groomer.png) 

Run FASTAQ groomer on the NA12878_CBW_chr1_R2.fastq file as well

File check:

![file3](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_3.png) 


#### 5) Trim the read and remove adapter sequence with Trimmomoatic
<a name="trim"></a>

** NGS: QC and manipulation/Trimmomatic   
 - Input FASTQ file: FASTQ groomer results for xxx_R1.fastq and xxx_R2.fastq files   
 - Perform initial ILLUMINACLIP step :Yes   
 -- Adapter sequence to use : TruSeq3 (paired-end, for MiSeq and HiSeq)   
 -- Maximum mismatch count which will still allow a full match to be performed : 2   
 -- How accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment: 30   
 -- How accurate the match between any adapter etc. sequence must be against a read: 15   
 - Select Trimmomatic operation to perform : insert Trimmomatic option   
 -- TRAILING: 20   
 -- MINLEN: 32   

![trim1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_Trim_1.png) 
![trim2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_Trim_2.png) 

This creates 4 files, for paired and unpaired reads for both R1 and R2

Use the paired results for the alignment

File check:

![file4](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_4.png) 


#### 6) Align the reads with BWA-MEM
<a name="align"></a>

** NGS: Mapping/Map with BWA-MEM 

 - Will you select a reference genome from your history or use a built-in index? : Use a genome from history or build index   
 - Use the following dataset as the reference sequence: hg19_chr1.fa   
 - Select first set of reads: Trimmomactic on FASTQ groomer (R1 paired)   
 - Select second set of reads: Trimmomactic on FASTQ groomer (R2 paired)   
 - Set the read groups information : Set read group SAM/BAM specification   
   with the following info  
 -- Read group identifier (ID): NA12878  
 -- Read group sample name (SM): NA12878   
 -- Platform/technology used to produce the reads (PL): ILLUMINA   
 -- Library name (LB):NA12878   
 -- Sequencing center that produced the read (CN): Broad Institute   
 -- Platform unit (PU): runNA12878_1      
 - Select Analysis node: 3. Full list of option   
 - Set input/output options : Set  
 - Mark shorter split hits of a chimeric alignment in the FLAG field as 'secondary alignment' instead of 'supplementary alignment': Yes   -M; For Picard<1.96 compatibility  

![align1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_align_1.png) 
![align2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_align_2.png) 

...

![align3](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_align_3.png) 

...

![align4](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_align_4.png) 

*** This step takes some time to run

File check:

![file5](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_5.png) 


#### 7) Sort the sam/bam 
<a name="sort"></a>

** NGS: Picard/SortSam sort SAM/BAM dataset  

  - Select SAM/BAM dataset or dataset collection: map with BWA-MEM file
  
  - Sort order: coordinate 
  
  - Select validation stringency: Silent   

![sort](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_sort.png) 

File check:

![file6](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_6.png) 


#### 8) Convert bam to sam file (optional)
<a name="convertBam"></a>

We will convert the bam file to a sam file to be able to look at it   
You can run this at different point of your analysis

** NGS: SAMtools/BAM-to-SAM

  - BAM File to Convert: sorted Bam file

![convertBam](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_convert_bam.png) 

File check:

![file7](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_7.png) 


#### 9) Indel realignment
<a name="indelRealign"></a>

** NGS GATK Tools/RealignerTargetCreator  

  - Choose the source for the reference list: History
  - Bam file: sorted Bam file
  - Using reference file: hg19_chr1.fa   
  - Basic or advance GATK option: Advanced   
  -- Operate on Genomic intervals: Genomic intervals : created interval on chr1
 
![realignertarget](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_realigner_target_1.png) 

...

![realignertarget](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_realigner_target_2.png) 

** NGS GATK Tools/IndelRealigner  

  - Choose the source for the reference list: History
  - Bam file: sorted Bam file
  - Using reference file: hg19_chr1.fa   
  - Restrict realignment to provided intervals: Realigner Target create results   

![indel](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_indel.png) 

File check:

![file8](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_8.png) 


#### 10) FixMates
<a name="fixmates"></a>

**NGS: Picard/FixMateInformation 

  - Select SAM/BAM dataset or dataset collection: Indel Realigner result bam
  - Select validation stringency: Silent   

![fixemate](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_fixemates.png) 

File check:

![file9](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_9.png) 


#### 11) Mark duplicates
<a name="markdup"></a>

**NGS: Picard/MarkDuplicates  

  - Select SAM/BAM dataset or dataset collection: FixMateInformation result   
  - Select validation stringency: Silent   

![markdup](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_markdup_1.png) 

...

![markdup2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_markdup_2.png) 


File check:

![file10](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_10.png) 


Have a look at the "Markduplicate metrics"


#### 12) Base Recalibration
<a name="recalibration"></a>

** NGS GATK Tools/BaseRecalibrator is not available!  

So we can use "Count covariates" and "Table recalibration". These two steps are the equivalent of BaseRecalibrator which is present in a newer version of GATK   

** NGS: GATK Tools/Count Covariates on BAM files

  - Choose the source for the reference list: History
  - Bam file: marked duplicates file   
  - Using reference genome: hg19_chr1.fa   
  - Covariates to be used in the recalibration: Select all, then unselect "Tile covariate"

![countCov](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_countcov.png) 


** NGS: GATK Tools/Table Recalibration on BAM files  

  - Covariates table recalibration file: Count Covariate 
  - Choose the source for the reference list: History
  - Bam File: Marked duplicates file  
  - Using reference genome:hg19_chr1.fa   

![table](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_table.png) 

File check:

![file11](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_11.png) 



#### 13) Extract Metrics
<a name="extracmetrics"></a>

** NGS GATK Tools/Depth of Coverage on BAM files   

 - Choose the source for the reference list: History
 - Bam file: Table recalibration result bam file
 - Using reference genome: hg19_chr1.fa   
 - Partition type for depth of coverage: select sample and readgroup
 - Summary coverage threshold   
  -- for summary file outputs, report the % of bases covered to >= this number: 10, 25, 50 and 100   **(insert 4 thresholds)**
 - Basic or Advanced GATK options: Advanced   
  -- Operate on Genomic intervals: Genomic intervals : created interval on chr1   
 - Basic or Advanced Analysis options: Advanced   
  -- "Omit the output of the depth of coverage at each base": Yes      

![cov1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_cov_1bis.png) 

...

![cov2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_cov_2.png) 

...

![cov3](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_cov_3.png) 

...

![cov5](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_cov_5.png) 

...

![cov4](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_cov_4.png) 


![file12a](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_12a.png) 

** NGS: Picard/Collect Alignment Summary Metrics writes a file containing summary alignment metrics 

 - SAM/BAM  dataset: Table recalibration result bam file   
 - Choose the source for the reference list: History   
 - Using reference genome: hg19_chr1.fa   
 - The level(s) at which to accumulate metrics set to: Read group, Sample      
 - Select validation stringency: Silent   

![sumA1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SumAlign_1.png) 

...

![sumA2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SumAlign_2.png) 

File check:

![file12b](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_12b.png) 

View "Collect Alignment Summary metrics"

** NGS: Picard/CollectInsertSizeMetrics plots distribution of insert sizes  

 - SAM/BAM  dataset: Table recalibration result bam file      
 - Choose the source for the reference list: History   
 - Using reference genome: hg19_chr1.fa   
 - The level(s) at which to accumulate metrics set to: Read group
 - Select validation stringency: Silent   

![sumI1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SumInsert_1.png) 

...

![sumI2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SumInsert_2bis.png) 

File check:

![file12c](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_12c.png) 

View "collectInsertSize Metrics" and pdf file



#### Variant calling and annotation from Module 5

#### All the following steps are detailed in the Module 5 practical

To continue you can use the aligned, sorted, marked duplicates and quality recalibrated files that you just created or download the one you used in Module 5 from the server (http://cbw##.dyndns.info/module5/NA12878.bwa.sort.rmdup.realign.bam, ## being your student id).


#### 14) Call SNPs   
<a name="callsnp"></a>

** NGS GATK Tools/Unified Genotyper SNP and indel caller  

 - Choose the source for the reference list: History   
 - BAM file: Table recalibration result bam file   
 - Using reference genome: hg19_chr1.fa         
 - The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold): 10   
 - Basic or Advanced GATK options: Advanced   
 -- Operate on Genomic intervals: Genomic intervals : created interval on chr1    

![uni1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_Unifier_1.png) 

...

![uni2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_Unifier_2.png) 

File check:

![file13](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_13.png) 


Have a look at the vcf file   


#### 15) Filter the variants  
<a name="filtervariant"></a>

Typically variant callers will only perform a minimal amount of filtering when presenting variant calls. In the case of GATK, we are actively removing any variant with score less than 10. Any variant with a score less than 30 is labeled with the filter “LowQual”.   

To perform more rigorous filtering, another program must be used. In our case, we will use the VariantFiltration tool in GATK.

NOTE: The best practice when using GATK is to use the VariantRecalibrator. In our data set, we had too few variants to accurately use the variant recalibrator and therefore we used the VariantFiltration tool instead.

** NGS GATK Tools/Variant Filtration on VCF files   

 - Choose the source for the reference list: History
 - Variant file to annotate: Unified genotyper results file 
 - Using reference genome: hg19_chr1.fa         
 - Variant filter   
 
set the 3 following filters  

 -- filter Expression: QD < 2.0 Filter, name: QDFilter   
 -- filter Expression: FS > 200.0 Filter, name: FSFilter    
 -- filter Expression: MQ < 40.0 Filter, name: MQFilter   

![varFil1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_VarFil_1.png) 

...

![varFil2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_VarFil_2.png) 

File check:

![file14](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_14.png) 


You can look at the output vcf file that contains "filter" annotation


#### 16) Annotate the variants
<a name="annotate"></a>

with snpEff   

** NGS: Variant Analysis/SnpEff Variant effect and annotation   

 - Sequence changes: Variant Filtration result   
 - Filter output  
  --select "Do not show INTERGENIC changes"   

![SnpEff1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SnpEff_1.png) 

...

![SnpEff2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_SnpEff_2.png) 

File check:

![file15b](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_15b.png) 

You can download the SnpEff stat html page and vcf file

with GATK VariantAnnotator 

** NGS GATK Tools/VariantAnnotator   

 - Choose the source for the reference list: History   
 - Variant file to annotate: SnpEff vcf file   
 - Using reference file: hg19_chr1.fa   
 - Provide a dbSNP reference-ordered data file: set dbSNP   
  -- ROD file: dbSNP_135_chr1.vcf.gz   
 - Basic or Advanced GATK options: Advanced   
 -- Operate on Genomic intervals: Genomic intervals : created interval on chr1    

![Anno1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_anno_1.png) 

...

![Anno2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_anno_2.png) 

...

![Anno4](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_anno_4.png) 

...

![Anno3](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_anno_3.png) 


File check:

![file15a](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_file_15a.png) 


Look at or download your filtered and annotated variant vcf files



#### The "History options" menu
<a name="menu"></a>

You can  

 - See saved histories
 - Extract workflow  
 - Save your workflow   
 - Share your workflow   
 - ...   

![history](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_history.png) 


#### Workflow
<a name="workflow"></a>

Use the workflow tab   
 - Edit your worklow   
 - Rename input and output files   
 - Run your workflow on other samples   
 ...

![workflow](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_workflow_2.png) 


For information (we will not do it as part of the practical), you can run a tool and/or a given workflow on several samples at once. You will need to "Build List of Dataset Pairs" then run the tools on "Dataset collection"

![s1](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_severalSample_1.png) 

 - Select your input samples 

![s2](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_severalSample_2.png) 

 - Create your pairs 
 - Run a tool on your dataset collection

![s3](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_severalSample_3.png) 

Next, you could...   
 - Create a new history, upload your input samples and run you workflow on them
 - Continue using the tools on Galaxy
 - Continue exploring/analysing/visualizing your data and results!

There are lots of tutorial, mainling list, videos for help  (such as https://vimeo.com/galaxyproject)   
**We hope that you enjoyed using Galaxy!**

