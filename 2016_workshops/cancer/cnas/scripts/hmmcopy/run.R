# Script to run HMMcopy which estimates copy number from genome sequencing data
# 
# Author: Andrew Roth
###############################################################################

## Setup path and sample information

# Set the directory where output will be placed
output.dir <- file.path(normalizePath('./'), 'results', 'hmmcopy')

# Name of sample we are analysing
sample.id <- 'HCC1143'

# Set the directory where the plots will be placed
plot.dir <- file.path(output.dir, 'plots')

# Set the directory where the results i.e. segment files will be placed
results.dir <- file.path(output.dir, 'results')

## Load the required utils

# Set working directory to where script is being run
setwd(dirname(sys.frame(1)$ofile))

# Load functions from helper scripts
source('../utils.R', chdir=TRUE)
source('plot.R')
source('utils.R')

# Load HMMcopy library
library(HMMcopy)

# Make the plot and results directories.
safe_make_dir(plot.dir)
safe_make_dir(results.dir)

#### Step 1 - Setup Input Files ####

# Path to file with GC content for each bin in WIG format 
gc.file <- '/home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.gc.wig'

# Path to file with mappability scores for each bin in WIG format.
map.file <- '/home/ubuntu/CourseData/CG_data/Module5/genome/hg19.21.map.wig'

# Path to file with normal read count data for each bin
normal.reads.file <- '/home/ubuntu/CourseData/CG_data/TCGA/HCC1143/G15511.HCC1143.1.chr21.wig'

# Path to file with tumour read count data for each bin
tumour.reads.file <- '/home/ubuntu/CourseData/CG_data/TCGA/HCC1143/G15511.HCC1143_BL.1.chr21.wig'

#### Step 2 - Read data into RangedData objects ####

# Read RangedData object for the normal data
normal.data <- wigsToRangedData(normal.reads.file, gc.file, map.file)

# Read RangedData object for the tumour data
tumour.data <- wigsToRangedData(tumour.reads.file, gc.file, map.file)

#### Step 3 - Correct read counts for GC content and mappability bias ####

# Correct read counts for the normal data
normal.corrected.data <- correctReadcount(normal.data)

# Correct read counts for the tumour data
tumour.corrected.data <- correctReadcount(tumour.data)

#### Step4 - Visualise the biases ####

# Path to file where the biases for the normal data will be plotted
normal.biases.file <- file.path(plot.dir, 'normal.biases.pdf') 

# Path to file where the biases for the tumour data will be plotted
tumour.biases.file <- file.path(plot.dir, 'tumour.biases.pdf')

# Plot the normal data biases
plotBiases(normal.corrected.data, normal.biases.file)

# Plot the tumour data biases
plotBiases(tumour.corrected.data, tumour.biases.file)

#### Step5 - Visualise the corrected read counts ####

# Path to file where the uncorrected and corrected read counts for the normal 
# data will be plotted
normal.readcount.file <- file.path(plot.dir, 'normal.readcount.pdf')

# Path to file where the uncorrected and corrected read counts for the tumour 
# data will be plotted
tumour.readcount.file <- file.path(plot.dir, 'tumour.readcount.pdf')

# Plot the read counts for the normal data
plotReadCounts(normal.corrected.data, normal.readcount.file)

# Plot the read counts for the normal data
plotReadCounts(tumour.corrected.data, tumour.readcount.file)

#### Step6 - Segment the tumour data ####

# Compute segments for the tumour data
tumour.segments <- HMMsegment(tumour.corrected.data)

# Inspect the first few segments
head(tumour.segments)

## First save the position (bin) level information ##

# Path to file where data for each position (bin) will be written
tumour.position.results.file <- file.path(results.dir, 'tumour.cna.tsv')

writeBinLevelResults(sample.id, tumour.corrected.data, tumour.segments,
		tumour.position.results.file)

## Next we save the segment level information in a format IGV can read ##

# Path where segment file will be written in IGV compatible format.
igv.segment.file <- file.path(results.dir, 'tumour.igv.seg')

# Save the results to file
writeSegmentResults(sample.id, tumour.corrected.data, tumour.segments,
		igv.segment.file, igv.compatible=TRUE)


## Finally we save the segment level information in a format APOLLOH can read ##			

# Path where segment file will be written in APOLLOH compatible format.
apolloh.segment.file <- file.path(results.dir, 'tumour.apolloh.seg')

# Save results to file
writeSegmentResults(sample.id, tumour.corrected.data, tumour.segments,
		apolloh.segment.file, igv.compatible=FALSE)
