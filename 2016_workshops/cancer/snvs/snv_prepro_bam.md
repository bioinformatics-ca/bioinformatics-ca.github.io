---
layout: post3
title: Lab Module 6 - Pre-processing Bams
header1: Bioinformatics for Cancer Genomics 2016
header2: Lab Module 6 - Pre-processing Bams
---

Often it is important to pre-process the bams to remove some poor quality reads to help increase the sensitivity of variant calling. Samtools provides functionality for this. For instance, you can remove PCR duplicates:

~~~bash
samtools rmdup in.bam out.bam
~~~

Or filter out reads that are not primary alignments (i.e align to multiple locations):

~~~bash
samtools filter -F 256 in.bam > out.bam
~~~

Or filter out reads that are not mapped:

~~~bash
samtools filter -F 4 in.bam > out.bam
~~~

The 4 and 256 flag is a bitwise flag that can be explained from this website.
