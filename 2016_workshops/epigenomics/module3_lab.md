---
layout: post2
permalink: /epigenomic_data_analysis_module3_lab_2016/
title: Epigenomic Data Analysis 2016 Student Page
header1: Epigenomic Data Analysis 2016
header2: Module 3 Lab
image: CBW_Epigenome-data_icon.jpg
---

# Module 3: Introduction to WGBS and Analysis 

## Important notes:
* Please refer to the following guide for instructions on how to connect to Guillimin and submit jobs: [using_the_guillimin_hpc.md](using_the_guillimin_hpc.md)
* The instructions in this tutorial will suppose you are in a Linux/Max environment. The equivalent tools in Windows are provided in the [Guillimin documentation](using_the_guillimin_hpc.md).
* The user **class99** is provided here as an example. You should replace it by the username that was assigned to you at the beginning of the workshop.


## Introduction

### Description of the lab
This module will cover the basics of Bisulfite-sequencing data analysis including data visualization in IGV.

### Local software that we will use
* ssh
* IGV


## Tutorial

### Getting started

#####  Connect to the Guillimin HPC
```
ssh class99@guillimin.clumeq.ca
```

You will be in your home folder. At this step, before continuing, please make sure that you followed the instructions in the section **"The first time you log in"** of the [Guillimin guide](using_the_guillimin_hpc.md). If you don't, compute jobs will not execute normally.

##### Prepare directory for module 4
```
rm -rf ~/module4
mkdir -p ~/module4
cd ~/module4
```

##### Copy data for module 4
```
cp /gs/project/mugqic/bioinformatics.ca/epigenomics/wgb-seq/data/* .
```

##### Check the files
By typing ```ls``` you should see something similar to this
```
[class99@lg-1r14-n04 module4]$ ls
fat.1.fastq  iPSC_1.1.fastq  iPSC_2.1.fastq  kidney.1.fastq
fat.2.fastq  iPSC_1.2.fastq  iPSC_2.2.fastq  kidney.2.fastq
```
*What do the ".1" and ".2" in the file names mean?*

### Map using bismark
We will now process and map the reads using Bismark.
```
echo 'module load mugqic/bismark/0.16.1 ; \
bismark --bowtie2 -n 1 /gs/project/mugqic/bioinformatics.ca/epigenomics/wgb-seq/genome/ \
-1 iPSC_1.1.fastq -2 iPSC_1.2.fastq' \
|  qsub -l nodes=1:ppn=1 -d .
```

##### Check files
```
[class99@lg-1r14-n04 module4]$ ls
fat.1.fastq	iPSC_1.1.fastq_C_to_T.fastq  iPSC_2.1.fastq  kidney.2.fastq
fat.2.fastq	iPSC_1.2.fastq		     iPSC_2.2.fastq  STDIN.e60334723
iPSC_1.1.fastq	iPSC_1.2.fastq_G_to_A.fastq  kidney.1.fastq  STDIN.o60334723
```
*Is this what you expected?*

##### Check for error message
```
less STDIN.e60334723
```
Where you replace the file name by your specific error file.

##### Map (again) using bismark
```
echo 'module load mugqic/bismark/0.16.1 ; module load mugqic/bowtie2/2.2.4 ; \
bismark --bowtie2 -n 1 /gs/project/mugqic/bioinformatics.ca/epigenomics/wgb-seq/genome/ \
-1 iPSC_1.1.fastq -2 iPSC_1.2.fastq' \
|  qsub -l nodes=1:ppn=1 -d .
```
