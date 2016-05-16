---
layout: post3
permalink: /bioinformatics_for_cancer_genomics_2016/rearrangement
title: Bioinformatics for Cancer Genomics 2016 Genome Rearrangement Tutorial
header1: Bioinformatics for Cancer Genomics 2016
header2: Genome Rearrangement Tutorial
image: CBW_cancerDNA_icon-16.jpg
---

# Introduction

In this tutorial we will use [lumpy-sv](https://github.com/arq5x/lumpy-sv) to perform the structural rearrangement analysis
using the BAM file from step 1.


# Genome Rearrangment Tutorial

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
    -pe id:sample,bam_file:sample.discordants.bam,histo_file:sample.lib1.histo,mean:330,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:sample,bam_file:sample.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > sample.vcf
```

The output is in VCF format: https://samtools.github.io/hts-specs/VCFv4.2.pdf

cat sample.vcf

Lumpy found about 260 SV events. The VCF file includes useful information about the variant calls.
The SVTYPE attribute tells us the structural variation type. Most of the events are deletions (SVTYPE=DEL).
The PE and SR tags in each line tell us how many paired-end and split-reads support each call. The calls with many
 supporting reads are usually more likely to be true structural variants.

Let's filter the list of events by only looking at those that are supported by 15 reads:

```
~/CourseData/CG_data/Module3/scripts/filter-by-support --min-support 15 sample.vcf
```

Only 44 events are supported by at least 15 reads. Let's now take a closer look in IGV.

## Viewing rearrangements in IGV and classifying events

# If you do not have IGV open, follow the instructions in lab 1.

# Navigate to the genomic location 20:32,407,585-32,416,459
# This region contains a deletion that lumpy-sv found.
# Right-click on an alignment and select 'view as pairs'
# If you hover or click on one of the reads that is colored red
# you can view its insert size. You will notice many of the pairs
# have an insert size of around 3000bp. These 'stretched' pairs
# support the deletion. The coverage track in IGV shows
# the deleted region has lower coverage, which also supports the call.

# Navigate to 20:33,113,802-33,118,238 to view a second deletion.

# Navigate to the location 20:25,967,655-25,976,529
# This region also shows 'stretched' pairs but there is no drop in coverage.
# If you examine the orientation of the pairs you will notice that they are both
# on the same strand (+/+ or -/-). Remember, we expect the pairs to be on opposite
# strands (+/- or -/+). This means the event is probably an inversion.

# Did lumpy find this call? Lets search the output for an event around this location
# We'll use awk to find all events on chromosome 29, between position 25,900,000 and 26,000,000

cat sample.vcf | awk '$1 == 20 && $2 >= 25900000 && $2 < 26000000'

# Lumpy found 5 events in this region. 4 of them are small deletions supported by less than 10 paired-end reads
# and no split reads. Lumpy also found a large inversion supported by 78 paired-end reads. This is the event we saw in IGV.
# It was not supported by split reads but this is not unusual for inversions. The paired-end evidence is high so
# we have confidence this is a true inversion.

# View the event at 20:26,318,238-26,318,791
# This region is difficult. Notice that the coverage significantly
# increases in the region and there are many reads mapped with mismatches
# and low mapping quality.
# Many of the pairs of the reads map to different chromosomes. It is difficult
# to tell what is happening in this case.

# For the remainder of the lab view other events and try to determine if they are real.

#
# End of lab 2
#
