---
layout: post2
permalink: /informatics_for_high-throughput_data_sequencing_2016_module7_lab/
title: Informatics for High-Throughput Sequencing Data 2016 Module 7 lab
header1: Informatics for High-Throughput Sequencing Data 2016
header2: Module 7 lab
image: CBW_High-throughput_icon.jpg
---
# CBW HT-seq Module 7 - Galaxy lab
This lab was created by Florence Cavalli.

## Table of contents
1. [Introduction](#introduction)
2. [Align DNA with BWA-MEM](#align)




## Introduction
<a name="introduction"></a>
This workshop will show you how to use Galaxy and performed the anaylis you ran in Module 2 and Module 5 with Galaxy.
We are therefore using the dataset from NA12878 patient (reads extracted from region chr1:17704860-18004860)

We will use the galaxy website: https://usegalaxy.org/
Start by loging on the website with your username

The paramaters that need to be set or changed are specified in the text below (and in the images), keep the rest of the parameters with their default values

** indicate the tool to use that you need to select from the left column
![tool](http://bioinformatics-ca.github.io/2016_workshops/ht-seq/img/Galaxy_tools.png) 

#### 1) Load the data 

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
#### 2) Convert the data to use Trimmomatic

** NGS: QC and manipulation/FASTAQ groomer

#### 3) Trim the read and remove adapter sequence with Trimmomoatic
** NGS: QC and manipulation/Trimmomatic 

- On FASTAQ groomer data 
- ILLUMINACLIP step
- TruSeq3 (paired-end, for MiSeq and HiSeq)
- 2
- 30
- 15
- TRAILING:20 
- MINLEN:32

It creates 4 files two for R1 and two for R2 (paired and unpaired?)
Use paired results for the alignment


#### 4) Align the reads with BWA-mem

** NGS: Mapping/Map with BWA-MEM 

- Set the Read group -RG option
- Set Read group information?
- Set read group SAM/BAM specification
with the following info
ID:NA12878
SM:NA12878
LB:NA12878
PU:runNA12878_1
CN:Broad Institute
PL:ILLUMINA

-Select Analysis node
 Full list of option
 Yes to -M option (for compatibility with Picard)


*** This step takes some time to run

#### 5) Sort the sam/bam 

** NGS: Picard/SortSam sort SAM/BAM dataset
- On aligned bam file
- Sort order: coordinate
- Select validation stringency: Silent

#### 6) Create single interval  

** Text manipulation/Create Single Interval
- Chr1:17704860-18004860

#### 7) Indel realignment

** NGS GATK Tools/RealignerTargetCreator 
- Bam file : Bam sorted in coordinate order
- Using reference file hg19_chr1.fa
- Basic or advance GATK option
 Advanced
 Operate on Genomic intervals -> Set -L parameter with your "Create single interval" data
 Keep other parameters as default


** NGS GATK Tools/IndelRealigner 
- Choose the source for the reference list: History
- Bam file : Bam sorted in coordinate order
- Using reference file hg19_chr1.fa

- Restrict realignment to provided intervals: Realigner Target create results


#### 8) FixMates

**NGS: Picard/FixMateInformation 
- Select SAM/BAM dataset or dataset collection ->Indel Realigner result
- Select validation stringency -> Silent
other default parameter

#### 9) Mark duplicates

**NGS: Picard/MarkDuplicates 
- Select SAM/BAM dataset or dataset collection -> FixMateInformation result
- Select validation stringency -> Silent

You can look at the Markduplicate metrics

#### 10) Recalibration
** NGS GATK Tools/BaseRecalibrator is not available!
So we can use "Count covariates" and "Table recalibration". These two step are teh equivalent of BaseRecalibrator which is present in a newer version of GATK

** NGS: GATK Tools/Count Covariates on BAM files
- On marked duplicat file
- Reference genome: hg19_Chr1.fa
- Select all the unselect Tile covariate


** NGS: GATK Tools/Table Recalibration on BAM files
- Covariates table recalibration file: Count Covariate
- Bam File: Mark duplicate
- Reference genome: hg19_Chr1.fa


#### 11) Extract Metrics
** NGS GATK Tools/Depth of Coverage on BAM files
- Summary coverage threshold
- insert 4 threshold at 10, 25, 50 and 100
- Basic or Advanced GATK options
  Advanced
Operate on Genomic intervals
-- L "create single interval

-Basic or Advanced Analysis options
 Advanced
"Omit the output of the depth of coverage at each base" set to Yes


** NGS: Picard/CollectInsertSizeMetrics plots distribution of insert sizes
- SAM/BAM  dataset -> Table recalibration BAM file
- Reference genome: hg19_Chr1.fa
- reference 
- The level(s) at which to accumulate metrics set to -> Library
- Select validation stringency -> silent


View collectInsertSize Metrics and pdf file


** NGS: Picard/Collect Alignment Summary Metrics writes a file containing summary alignment metrics
- SAM/BAM  dataset -> Table recalibration BAM file
- Reference genome: hg19_Chr1.fa
- The level(s) at which to accumulate metrics set to -> Library
- Select validation stringency -> silent

View Collect Alignment Summary metrics


Variant calling an annotation from Module 5

Use your the aligned, sorted  with removed duplicates files that you just created (or download the one you used in Module 5  from the server (NA12878.bwa.sort.rmdup.realign.bam))

#### 12) Call SNPs in NA12878 using the sorted, realigned removed duplicates bam files
** NGS GATK Tools/Unified Genotyper SNP and indel caller
- BAM file -> marked duplicates
- reference genome -> hg19_chr1.fa
- The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold): 10

- Basic or Advanced GATK options
 Advanced
 Operate on Genomic intervals
-- L "create single interval

Have a look at the vcf file


#### 13) Filter the variants

Typically variant callers will only perform a minimal amount of filtering when presenting variant calls. In the case of GATK, we are actively removing any variant with score less than 10. Any variant with a score less than 30 is labeled with the filter “LowQual”.

To perform more rigorous filtering, another program must be used. In our case, we will use the VariantFiltration tool in GATK.

NOTE: The best practice when using GATK is to use the VariantRecalibrator. In our data set, we had too few variants to accurately use the variant recalibrator and therefore we used the VariantFiltration tool instead.


** NGS GATK Tools/Variant Filtration on VCF files
- BAM file -> marked duplicates
- Reference genome -> hg19_chr1.fa

- Variant filter
set the 3 following filters
filter Expression:QD < 2.0 Filter name:QDFilter
filter Expression:FS > 200.0" Filter name:FSFilter 
filter Expression:MQ < 40.0 Filter name:MQFilter


You can look at the output vcf file that contains some filter annotation

#### 14) Annotate the vcf file
with snpEff
** NGS: Variant Analysis/SnpEff Variant effect and annotation
- Sequence changes: Variant Filtration result
- Both input and output as vcf file (default)
- Genome source: GRCh37.74:hg19
- Filter out
  select "Do not show INTERGENIC changes"

Look at the output (vcf file and stats)


with GATK VariantAnnotator 

** NGS GATK Tools/VariantAnnotator
- Variant file to annotate: SnpEff file
- reference genome: hg19_chr1.fa
-Provide a dbSNP reference-ordered data file
 set dbSNP: dbSNP_135_chr1.vcf.gz
- Basic or Advanced GATK options
Advanced
Operate on Genomic intervals
-- L -> create single interval"

Look at or download your filtered and annotated variant calls as vcf files


#### Look at the "History options" menu
- Extract workflow
- Save your workflow
- Share your workflow
...

- continue using the tools on galaxy…
- continue exploring/analysing/visualizing your data and results..


There are lots of tutorial, mainling list, videos for help  (such as https://vimeo.com/galaxyproject)/


