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

bwa mem is the leading algorithm for mapping short Illumina reads to a reference genome. You can read about how it works [here](http://arxiv.org/abs/1303.3997). We'll now run bwa mem to map HCC1395 reads to the human reference genome. 

In the following command we provide bwa with the location of the reference genome - in this exercise we use the human reference genome prepared by the [1000 Genomes project](http://www.1000genomes.org/category/reference/) - and a FASTQ file containing the paired end reads for the tumour sample. In this case the paired end reads are in *interleaved* format where the two halves of a pair are in consecutive FASTQ records. The output file is in the SAM format. When mapping a whole-genome sequencing run it will take many hours to run - in this tutorial we only use a subset of the reads so this step doesn't take very long.

```
bwa mem -t 4 -p ~/CourseData/CG_data/Module3/human_g1k_v37.fasta ~/CourseData/CG_data/Module3/reads.tumour.fastq > tumour.sam
```


## Exploring the alignments

SAM is a plain-text format that can be viewed from the command line. You can use the head command to look at the alignments:

```
head -100 tumour.sam
```

You will see the SAM header containing metadata followed by a few alignments. You can refer to the slides from the lecture to determine the meaning of each field.

In this SAM file, the reads are ordered by their position in the original FASTQ file. Most programs want to work with the alignments ordered by their position on the reference genome. We'll use [samtools](https://github.com/samtools/samtools) to sort the alignment file. To do this, we need to first convert the SAM (text) to the BAM (binary) format. We use the `samtools view -Sb` command to do this, and pipe the output directly into samtools sort.

```
samtools view -Sb tumour.sam | samtools sort -o tumour.sorted.bam
```

The samtools view command can also be used to convert BAM to SAM:

```
samtools view tumour.sorted.bam | head -100
```

samtools also provides functions to view the alignments for a particular region of the reference genome. To do this we first need to build an index of the BAM file, which allows samtools to quickly extract the alignments for a region without reading the entire BAM file:

```
samtools index tumour.sorted.bam
```

Now that the BAM is indexed, we can view the alignments for any region of the genome:

```
samtools view tumour.sorted.bam 9:14,196,000-14,197,000
```

As a tab-delimited file, the SAM format is easy to manipulate with common unix tools like grep, awk, cut, sort. For example, this command uses cut to extract just the reference coordinates from the alignments:

```
samtools view tumour.sorted.bam 9:14,196,000-14,197,000 | cut -f3-4
```


The `samtools flagstat` command gives us a summary of the alignments:

```
samtools flagstat tumour.sorted.bam
```

We can use `samtools idxstats` to count the number of reads mapped to each chromosome:

```
samtools idxstats tumour.sorted.bam
```

This command will just display the number of reads mapped to chromosome 9:

```
samtools idxstats tumour.sorted.bam | awk '$1 == "9"'
```

## Examining alignments

You can view all of the aligned reads for a particular reference base using mpileup.
The following command will show the read bases and their quality scores at a heterozygous SNP.
The "." and "," symbols indicate bases that match the reference. 
There are 38 reads that show a "G" base at this position. 
The individual's genotype at this position is likely A/G.

```
samtools mpileup -f ~/CourseData/CG_data/Module3/human_g1k_v37.fasta -r 9:14,196,087-14,196,087 tumour.sorted.bam
```

Load the data into IGV by performing the following:

```
   Open IGV and change the genome from hg19 to 'human_g1k_v37'
   Choose 'Load from URL' from the file menu
   Type: http://cbw#.dyndns.info/Module3/tumour.sorted.bam where # is your student ID
   Navigate to 9:14,196,087
```

Notice that the alignments have high mapping quality.

This is the end of the mapping tutorial. In the remaining time you can also map reads for the matched-normal sample for this cell line. The reads are in `~/CourseData/CG_data/Module3/reads.normal.fastq`.
