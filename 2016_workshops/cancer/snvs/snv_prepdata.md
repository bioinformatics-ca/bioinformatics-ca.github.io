---
layout: post3
title: Lab Module 6 - Data Preparation
header1: Bioinformatics for Cancer Genomics 2016
header2: Lab Module 6 - Data Preparation
---

## Introduction

This lab uses a number of publicly available datasets.  Below we give details about how to obtain these datasets.  We also include instructions for the time intensive processing has been applied to these datasets prior to the lab.  The instructions below will not be needed during the lab, they are included for reference purposes in case you wish to reproduce the lab in your own time.

## Preparing the Sequencing Data

We will restrict our analysis to a 1 Mb region (7Mb and 8Mb) within chromosome 17. These bam files

Use `samtools` to create bam files containing only alignments within this region.  The argument `17:7000000-8000000` specifies the region, and `-b` specifies the output is bam (default output is sam format).

~~~bash
mkdir data
samtools view -b HCC1395/exome/HCC1395_exome_tumour.bam 17:7000000-8000000 \
    > data/HCC1395_exome_tumour.17.7MB-8MB.bam
samtools view -b HCC1395/exome/HCC1395_exome_normal.bam 17:7000000-8000000 \
    > data/HCC1395_exome_normal.17.7MB-8MB.bam
~~~

Create an index for each bam file.

~~~bash
samtools index data/HCC1395_exome_tumour.17.7MB-8MB.bam
samtools index data/HCC1395_exome_normal.17.7MB-8MB.bam
~~~
