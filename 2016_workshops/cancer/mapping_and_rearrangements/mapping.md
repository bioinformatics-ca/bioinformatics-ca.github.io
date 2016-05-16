---
layout: post3
permalink: /bioinformatics_for_cancer_genomics_2016/mapping
title: Bioinformatics for Cancer Genomics 2016 Read Mapping Tutorial
header1: Bioinformatics for Cancer Genomics 2016
header2: Read Mapping Tutorial
image: CBW_cancerDNA_icon-16.jpg
---

# Introduction

In this tutorial we will use [bwa](https://github.com/lh3/bwa) to map reads to the human reference genome.

The data that we will use is from a breast cancer cell line, HCC1395. In the second part of this module we will use [lumpy](https://github.com/arq5x/lumpy-sv) to find genome rearrangements using the mapped reads. Module 4 will show you how to discover fusion genes using the same data set.

# Mapping Tutorial


## Preparation

First, let's create a directory to hold our results:

```
cd ~/workspace
mkdir Module3
cd Module3
```

## Mapping using bwa-mem

bwa mem is the leading mapping algorithm for Illumina short reads. You can read about how it works [here](http://arxiv.org/abs/1303.3997). We'll now run bwa mem to map HCC1395 reads to the human reference genome. We give bwa mem two FASTQ files, one for the first half of the pair and one for the second half of the pair. The output file is in the SAM format. When mapping a whole-genome sequencing run it will take many hours to run - in this tutorial we'll only use a subset of the reads that will map quickly.

```
bwa mem -t 4 ~/CourseData/CG_data/Module3/human_g1k_v37.fasta ~/CourseData/CG_data/Module3/HCC1395_sample_1.fastq ~/CourseData/CG_data/Module3/HCC1395_sample_2.fastq > sample.sam
```


## Exploring the alignments

SAM is a plain-text format that can be viewed from the command line. You can use the head command to look at your alignments:

```
head -100 sample.sam
```

You will see the SAM header containing metadata followed by a few alignments.

In this SAM file, the reads are ordered by their position in the original FASTQ file. Most programs want to work with the alignments ordered by their position on the reference genome. We'll use samtools to sort the alignment file. To do this, we need to first convert the SAM (text) to the BAM (binary) format. We use the `samtools view -Sb` command to do this, and pipe the output directly into samtools sort.

```
samtools view -Sb sample.sam | samtools sort -o sample.sorted.bam
```

The samtools view command can also be used to convert BAM to SAM

```
samtools view sample.sorted.bam | head -100
```

samtools also provides functions to request the alignments for particular regions of the reference genome. To do this we first need to build an index of the BAM file, which allows samtools to quickly extract the alignments for a region without reading the entire BAM file:

```
samtools index sample.sorted.bam
```

Now that the BAM is indexed, we can view the alignments for any region of the genome:

```
samtools view sample.sorted.bam 20:26,000,000-26,010,000
```

As a tab-delimited file, the SAM format is easy to manipulate with common unix tools like grep, awk, cut, sort. For example, this command uses cut to extract just the reference coordinates from the alignments:

```
samtools view sample.sorted.bam 20:26,000,000-26,010,000 | cut -f3-4
```


The `samtools flagstat` command gives us a summary of the alignments:

```
samtools flagstat sample.sorted.bam
```

We can use `samtools idxstats` to count the number of reads mapped to each chromosome:

```
samtools idxstats sample.sorted.bam
```

This command will just display the number of reads mapped to chromosome 20

```
samtools idxstats sample.sorted.bam | awk '$1 == "20"'
```

## Examining alignments

You can view all of the aligned reads for a particular reference base using mpileup.
The following command will show the read bases and their quality scores at a heterozygous SNP.
The "." and "," symbols indicate bases that match the reference. 
There are 18 reads that show a "G" base at this position. 
The individual's genotype at this position is likely A/G.

```
samtools mpileup -f ~/CourseData/CG_data/Module3/human_g1k_v37.fasta -r 20:32,001,292-32,001,292 sample.sorted.bam
```

Load the data into IGV by performing the following:

```
   Open IGV and change the genome from hg19 to 'human_g1k_v37'
   Choose 'Load from URL' from the file menu
   Type: http://cbw#.entrydns.org/module3/sample.sorted.bam where # is your student ID
   Navigate to 20:32,001,292
```

Notice that the alignments have high mapping quality.

Now we will view the alignments for a different position:

```
samtools mpileup -f ~/CourseData/CG_data/Module3/human_g1k_v37.fasta -r 20:25,997,273-25,997,273 sample.sorted.bam
```

In this case, 11 reads show a T base at this position.
Look at the alignments in IGV by navigating to 20:25,997,273. 
The reads colored white have mapping quality 0.
This means the alignment is ambiguous and should not be trusted.
It is unclear whether the T->C alignments are true SNPs.

This is the end of lab 1. You can use the remaining time to explore the alignments
and ask questions if you notice anything unusual or interesting.

