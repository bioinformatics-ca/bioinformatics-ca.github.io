---
layout: post2
permalink: /epigenomic_data_analysis_module1_lab_2016/
title: Epigenomic Data Analysis 2016 Student Page
header1: Epigenomic Data Analysis 2016
header2: Module 1 Lab
image: CBW_Epigenome-data_icon.jpg
---

# Module 1: Introduction to ChIP sequencing & analysis 

## Important notes:
* Please refer to the following guide for instructions on how to connect to Guillimin and submit jobs: [using_the_guillimin_hpc.md](using_the_guillimin_hpc.md)
* The instructions in this tutorial will suppose you are in a Linux/Max environment. The equivalent tools in Windows are provided in the [Guillimin documentation](using_the_guillimin_hpc.md).
* The user **class99** is provided here as an example. You should replace it by the username that was assigned to you at the beginning of the workshop.


## Introduction

### Description of the lab
This module will cover the basics of how to login to the cluster, launch jobs and perform a basic QC analysis of the data sets.

### Local software that we will use
* ssh
* Web browser to visualize FastQC output


## Tutorial

### Getting started

#####  Connect to the Guillimin HPC
```
ssh class99@guillimin.clumeq.ca
```

You will be in your home folder. At this step, before continuing, please make sure that you followed the instructions in the section **"The first time you log in"** of the [Guillimin guide](using_the_guillimin_hpc.md). If you don't, compute jobs will not execute normally.

##### Prepare directory for module 1
```
rm -rf ~/module1
mkdir -p ~/module1
cd ~/module1
```

### Assessing FASTQ file quality with FastQC

##### Copy locally the FASTQ file that we will need for our FastQC analysis
```
cp /gs/project/mugqic/bioinformatics.ca/epigenomics/chip-seq/H1/data/H3K27ac/H3K27ac.H1.fastq.gz .
```

##### Check files

At this point if you type ```ls``` should have something like:
```
[class99@lg-1r14-n04 module1]$ ls
H3K27ac.H1.fastq.gz
```

#####  Run the FastQC command on the scheduler
```
echo 'module load mugqic/fastqc/0.11.2 ; fastqc H3K27ac.H1.fastq.gz' | qsub -l nodes=1:ppn=1 -d .
```

#####  Check the status of the job
```
showq -uclass%%
```
Where you replace **%%** by your student number. It usually takes a few seconds/minutes for the job to appear depending on the load of the cluster.

##### Check files

At this point if you type ```ls``` should have something like
```
[class99@lg-1r14-n04 module1]$ ls
H3K27ac.H1_fastqc.html	H3K27ac.H1_fastqc.zip  H3K27ac.H1.fastq.gz  STDIN.e60293217  STDIN.o60293217
```

#####  Download the results to your local computer
```
scp class99@guillimin.clumeq.ca:/home/class99/module1/H3K27ac.H1_fastqc.html .
```

#####  Open the downloaded file in a web browser

Open the folder and then double-click the file 
```
open .
```

Or directly from the command line using a command such as
```
firefox H3K27ac.H1_fastqc.html
```
