---
layout: post3
permalink: /bioinformatics_for_cancer_genomics_2016/rearrangement
title: Bioinformatics for Cancer Genomics 2016 Genome Rearrangement Tutorial
header1: Bioinformatics for Cancer Genomics 2016
header2: Genome Rearrangement Tutorial
image: CBW_cancerDNA_icon-16.jpg
---

## Introduction

In this tutorial we will use [lumpy-sv](https://github.com/arq5x/lumpy-sv) to perform the structural rearrangement analysis
using the BAM file for the HCC1395 reads we mapped in step 1.

## Preparing Input

We prepare the input for lumpy by making BAM files that only contain discordant read pairs and split reads.

First, we will make a BAM containing discordant pairs. Remember, discordant read pairs are those that do not map in the expected orientation or are too close/too far apart.

```
samtools view -b -F 1294 tumour.sorted.bam > tumour.discordants.bam
```

Next, make a BAM containing split reads. Split reads are those that are mapped with a large insertion or deletion in the alignment.

```
samtools view -h tumour.sorted.bam | ~/CourseData/CG_data/Module3/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > tumour.splitters.bam
```

Determine the distribution of paired-end fragment sizes. This tells lumpy how far apart paired-end reads should be and helps it determine where structural variants are. This code will print the mean and standard deviation of the fragment distribution to the screen.

```
cat tumour.sam | tail -n+100000 | ~/CourseData/CG_data/Module3/scripts/pairend_distro.py -r 101 -X 4 -N 10000 -o sample.lib1.histo
```

## Running Lumpy

Now we are ready to run lumpy. The following command is long so we break it into multiple lines by using the "\" symbol. You should still be able to copy-and-paste it into your terminal. We are passing the discordant read BAM file in with the -pe flag, and the split read BAM file in with the -se flag. 

```
lumpy \
    -mw 4 \
    -tt 0 \
    -pe id:tumour,bam_file:tumour.discordants.bam,histo_file:sample.lib1.histo,mean:250,stdev:40,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:tumour,bam_file:tumour.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > tumour.vcf
```

The output is in [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

```
cat tumour.vcf
```

Lumpy found about 50 rearrangement events in the small sample of reads that we are using for the tutorial. The VCF file includes useful information about the variant calls.  The SVTYPE attribute tells us the structural variation type. Many of the events are deletions (SVTYPE=DEL).  The PE and SR tags in each line tell us how many paired-end and split-reads support each call. The calls with many
 supporting reads are usually more likely to be true structural variants.


## Viewing rearrangements in IGV and classifying events

If you do not have IGV open, follow the instructions in lab 1.

Navigate to the genomic location `9:14,205,626-14,206,175`.
This region contains a small deletion that lumpy-sv found.
You can right-click on an alignment and select 'view as pairs' to have IGV draw a line linking the two ends of a read.
You will notice many of the pairs have a larger-than-expected insert size - IGV colors these pairs red. 
These 'stretched' pairs support the deletion. Also, there are quite a few reads whose alignment ends right at the breakpoint, which also supports the deletion, as does the aburpt drop in coverage.
It is very likely this is a true deletion.

Navigate to `9:108,329,845-108,347,463` to view a second, larger, deletion.

Navigate to the location `6:89,554,173-89,554,839`. In this case there are a number of pairs colored blue that indicate the other half of the pair maps to a different chromosome. If you click on one of these pairs you will see the other half maps to chromosome 1. There are many such pairs, and they have high mapping quality, which suggests this event might be a true rearrangement. If you right-click on one of the colored pairs and select "Go to mate" IGV will jump to the corresponding region on chromosome 1. In the coverage track you will notice that the read depth changes at the breakpoint on both chromosome 6 and chromosome 1. This suggests the rearrangement might also involve a copy number abnormality.

If you load the BAM of the normal sample you will notice there is no copy number change and no paired reads indicating the rearrangement in the normal sample, suggesting this is a somatic event - a change that only occurs in the tumour genome.  At location `12:24,104,965-24,106,007` there is another example of a somatic genome rearrangement.

Naviate to the location `4:12,098,102-12,104,063` to view an example of a small inversion - note the orientation of the pairs as they are both on the `-` strand (`-/-` orientation) whereas we expect normal pairs to be `+/-`. At `3:80,227,206-80,230,293` there is an example of a much larger inversion.

Finally view the rearrangement at `6:46,608,429-46,609,095`, this event is a gene fusion - it will be discussed later in Module 4.

For the remainder of the lab view other events that lumpy found and try to determine what type of rearrangement they are.
