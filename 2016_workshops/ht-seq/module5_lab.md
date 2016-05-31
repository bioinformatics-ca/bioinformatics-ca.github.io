---
layout: post2
permalink: /informatics_for_high-throughput_data_sequencing_2016_module5_lab/
title: Informatics for High-Throughput Sequencing Data 2016 Module 5 lab
header1: Informatics for High-Throughput Sequencing Data 2016
header2: Module 5 lab
image: CBW_High-throughput_icon.jpg
---

This lab was created by Guillaume Bourque.

## Introduction

In one of the previous modules, we have aligned the reads from NA12878 (daughter) in a small region on chromosome 1. In this module, we will use the data in the BAM files to call variants and then we'll perform some annotation and filtration.

As a brief reminder, these are the IDs for each member in the trio:
<table><thead>
<tr>
<th>Mother</th>
<th>Father</th>
<th>Child</th>
</tr>
</thead><tbody>
<tr>
<td>NA12892</td>
<td>NA12891</td>
<td>NA12878</td>
</tr>
</tbody></table>

Our read were extracted from the following regions:
<table><thead>
<tr>
<th>Chromosome</th>
<th>Start</th>
<th>End</th>
</tr>
</thead><tbody>
<tr>
<td>chr1</td>
<td>17704860</td>
<td>18004860</td>
</tr>
</tbody></table>

### Preliminaries

#### Amazon node

Read these directions [add link] for information on how to log in to your assigned Amazon node. 

#### Work directory

Create a new directory that will store all of the files created in this lab.

<pre><code>rm -rf ~/workspace/module5
mkdir -p ~/workspace/module5
cd ~/workspace/module5
ln -s ~/CourseData/HT_data/Module5/* .
</code></pre> 

***Note:***
    The `ln -s` command adds symbolic links of all of the files contained in the (read-only) `~/CourseData/HT_data/Module5` directory.
    
#### Input files

Our starting data set consists of 100 bp paired-end Illumina reads from the child (NA12878) that have been aligned to hg19 during one of the previous modules (NA12878.bwa.sort.bam). We also have the same data after duplicate removal and realignment around indels NA12878.bwa.sort.rmdup.realign.bam.

If you type ls, you should have something like:

<pre><code>ubuntu@ip-10-182-231-187:~/workspace/module5$ ls
NA12878.bwa.sort.bam      NA12878.bwa.sort.rmdup.realign.bai  other_files
NA12878.bwa.sort.bam.bai  NA12878.bwa.sort.rmdup.realign.bam
</pre></code>
`other_files` is a directory that contains additional reference files and scripts that we will need during the module.

Do you know what are the `.bai` files?

If you're interested, you can also look at the fastQC reports for the original fastq files in the directory:

`~workspace/module5/other_files/fastQC_reports`

You can view these reports at http://cbw##.dyndns.info/module5/other_files/fastQC_reports/.

***Note:*** You need to replace ## by your student number. 

## Calling variants with GATK

If you recall from the previous module, we first mapped the reads to hg19 and then we removed duplicate reads and realigned the reads around the indels.

Let's call SNPs in NA12878 using both the original and the improved bam files:

<pre><code>java -Xmx2g -jar /usr/local/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -l INFO \
 -R other_files/hg19.fa -I NA12878.bwa.sort.bam -stand_call_conf 30 -stand_emit_conf 10 \
 -o NA12878.bwa.sort.bam.vcf -nt 4 -glm BOTH -L chr1:17704860-18004860
java -Xmx2g -jar /usr/local/GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper -l INFO \
-R other_files/hg19.fa -I NA12878.bwa.sort.rmdup.realign.bam -stand_call_conf 30 \
-stand_emit_conf 10 -o NA12878.bwa.sort.rmdup.realign.bam.vcf -nt 4 \
-glm BOTH -L chr1:17704860-18004860
</pre></code>

`-Xmx2g` instructs java to allow up 2 GB of RAM to be used for GATK.
 
`-l` INFO specifies the minimum level of logging. 
 
`-R` specifies which reference sequence to use. 
 
`-I` specifies the input BAM files. 

`-stand_call_conf` is the minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls (that aren't at 'trigger' sites). Only genotypes with confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 (this is the default). 

`-stand_emit_conf` is the minimum phred-scaled Qscore threshold to emit low confidence calls (that aren't at 'trigger' sites). Genotypes with confidence >= this but less than the calling threshold are emitted but marked as filtered. The default value is 30. 

`-L` indicates the reference region where SNP calling should take place 

`-nt` 4 specifies that GATK should use 4 processors. 

`-glm` BOTH instructs GATK that we want to use both the SNP and INDEL genotype likelihoods models 

#### File check

At this point, you should have the following result files 

<pre><code>ubuntu@ip-10-182-231-187:~/workspace/module5$ ls
NA12878.bwa.sort.bam      NA12878.bwa.sort.bam.vcf.idx        NA12878.bwa.sort.rmdup.realign.bam.vcf
NA12878.bwa.sort.bam.bai  NA12878.bwa.sort.rmdup.realign.bai  NA12878.bwa.sort.rmdup.realign.bam.vcf.idx
NA12878.bwa.sort.bam.vcf  NA12878.bwa.sort.rmdup.realign.bam  other_files
</pre></code>

## Investigating the SNP calls

Use less to take a look at the vcf files:

<pre><code>less NA12878.bwa.sort.bam.vcf
less NA12878.bwa.sort.rmdup.realign.bam.vcf
</pre></code>

vcf is a daunting format at first glance, but you can find some basic information about the format here or here. [add links]

How do you figure out what the genotype is for each variant?

Do we have any annotation information yet?

How many SNPs were found?

Did we find the same number of variants using the files before and after duplicate removal and realignment? 