# Module 1: Introduction to ChIP sequencing & analysis 

## Important notes:
* Please refer to the following guide for instructions on how to connect to Guillimin and submit jobs: [using_the_guillimin_hpc.md](using_the_guillimin_hpc.md)
* The instructions in this tutorial will suppose you are in a Linux/Max environment. The equivalent tools in Windows are provided in the documentation [here](using_the_guillimin_hpc.md).
* The user **class99** is provided here as an example. You should replace it by the username that was assigned to you at the beginning of the workshop.


## Introduction

### Description of the lab
What will be covered in this module

### Software that we will use
* ssh
* Web browser to visualize FastQC output


## Tutorial

### Assessing FASTQ file quality with FastQC
* Copy locally the FASTQ file that we will need for our FastQC analysis.
```
cp /gs/project/mugqic/bioinformatics.ca/epigenomics/chip-seq/H1/data/H3K27ac/H3K27ac.H1.fastq.gz .
```

* Run the FastQC command on the scheduler.
```
echo 'module load mugqic/fastqc/0.11.2 ; fastqc H3K27ac.H1.fastq' | qsub -l nodes=1:ppn=1 -d .
```

* Download the results to your local computer
```
scp class99@guillimin.clumeq.ca:/home/class99/H3K27ac.H1_fastqc.html .
```

* Open the downloaded file in a web browser, either by double-clicking the file, or from the command line using a command such as
```
firefox H3K27ac.H1_fastqc.html
```