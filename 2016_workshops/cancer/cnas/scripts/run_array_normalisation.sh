#!/usr/bin/env bash                                                                    
####### NORMALIZE AFFYMETRIX SNP6 DATA USING PENNCNV-AFFY #####################################
####### REFER TO http://www.openbioinformatics.org/penncnv/penncnv_tutorial_affy_gw6.html #####
###############################################################################################

####### YOU NEED TO DOWNLOAD AND INSTALL THE FOLLOWING ################################################
### 1) http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz #######################
### 2) http://www.openbioinformatics.org/penncnv/download/gw6.tar.gz ##################################
### 3) http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip ##
### 4) http://www.affymetrix.com/support/developer/powertools/index.affx ##############################
#######################################################################################################

### Set up variables ###
rawCelPath=$1  #input directory containing your batch of SNP6 CEL files

outDir=$2  #output directory for your results

# Path to Affymetrix Power Tools
aptDir=/usr/local/apt

# Path to Affymetrix GenomeWide SNP6 files distributed from PennCNV
gw6Dir=~/CourseData/CG_data/Module5/data/config/gw6

# Path Affymetrix GenomeWide SNP6 files from affymetrix
affyDir=~/CourseData/CG_data/Module5/data/config/CD_GenomeWideSNP_6_rev3

# Path to PennCNV installation
penncnvDir=/usr/local/penncnv

## The following files are all required by the normalisation software
## These should not need to be changed for analysis of SNP6 arrays
cdfFile=$affyDir/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.cdf

sketchFile=$gw6Dir/lib/hapmap.quant-norm.normalization-target.txt

clusterFile=$gw6Dir/lib/hapmap.genocluster

locFile=$gw6Dir/lib/affygw6.hg18.pfb

### Step 1: PROBESET SUMMARIZATION ####################################################################
### Must have these installed: ########################################################################
### 4) http://www.affymetrix.com/support/developer/powertools/index.affx ##############################
### 3) http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip ##
echo "Running Step 1 - probeset summarization"
$aptDir/bin/apt-probeset-summarize --cdf-file $cdfFile --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch $sketchFile --out-dir $outDir/apt --cel-files $rawCelPath --chip-type GenomeWideEx_6 --chip-type GenomeWideSNP_6

### Step 2: PROBESET SUMMARIZATION ####################################################################
### Must have these installed: ########################################################################
### 1) http://www.openbioinformatics.org/penncnv/download/penncnv.latest.tar.gz########################
### 3) http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip ##
echo "Running Step 2 - B-allele and log ratios"
$gw6Dir/bin/normalize_affy_geno_cluster.pl $clusterFile  $outDir/apt/quant-norm.pm-only.med-polish.expr.summary.txt -locfile $locFile -out $outDir/gw6.lrr_baf.txt

### Step 3: SPLIT RESULTS TO INDIVIDUAL FILES
### The last command puts the results for all arrays in the batch, into the same file. For subsequent 
### analysis it is useful to create a file for each array. This command creates the individual files.
echo "Running Step 3 - Splitting results into separate files."
$penncnvDir/kcolumn.pl $outDir/gw6.lrr_baf.txt split 2 -tab -head 3 -name --output $outDir/gw6
