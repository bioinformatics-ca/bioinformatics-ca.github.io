---
layout: post3
title: Bioinformatics for Cancer Genomics 2016 Gene Fusions Tutorial
header1: Bioinformatics for Cancer Genomics 2016
header2: Gene Fusions Tutorial
image: fusions.png
---

# General installation for CBW tutorial

## Environment setup

For this tutorial, we now assume the environment variable `$TUTORIAL_HOME`
has been set to an existing directory to which the user has write access.

All binaries used in this tutorial will be installed using [`conda`](https://www.continuum.io/downloads).  Modify the `PATH` environment variable to point to binaries in the anaconda installation.

~~~
export PATH="/home/ubuntu/CourseData/CG_data/Module4/anaconda/bin:$PATH"
~~~

 

## Installation

### Tutorial directory structure

Create a directory for the reference data.

~~~ bash
mkdir -p $TUTORIAL_HOME/refdata
~~~

### Install tutorial scripts

All content for this tutorial is located in a bitbucket repo at
`https://dranew@bitbucket.org/dranew/cbw_tutorial.git`.  Some of the data, scripts and config files from this repo will be used in the tutorial.  Clone the tutorial repo so we have a copy of the tutorial scripts in a known location.

~~~ bash
cd $TUTORIAL_HOME/
~~~

~~~ bash
git clone https://bitbucket.org/dranew/cbw_tutorial.git
~~~

### Install Anaconda

All binaries used in this tutorial will be installed using [`conda`](https://www.continuum.io/downloads).  Download and install anaconda with the prefix `$TUTORIAL_HOME/anaconda`.

Packages in `conda` are stored in channels.  Several additional channels hosting bioinformatics specific software will be required.  Add additional channels using `conda config`.

~~~
conda config --add channels r
conda config --add channels bioconda
conda config --add channels BioBuilds
conda config --add channels https://conda.anaconda.org/dranew
~~~

### Install samtools

The `samtools` package is the most widely used software for manipulating high-throughput sequence data stored in the 'bam' format.

Install `samtools` using `conda`.

~~~ bash
conda install samtools
~~~

### Install picard tools

Picard tools is a useful set of utilities for manipulating sequence data in bam/sam format.

Install `picard` using `conda`.

~~~ bash
conda install picard
~~~

### Install igv tools

The igvtools package provides utilities for preprocessing bam files for quicker viewing in IGV.

Install `igvtools` using `conda`.

~~~ bash
conda install igvtools
~~~

### Install bowtie and bowtie2

Install `bowtie` and `bowtie2` using `conda`.

~~~ bash
conda install bowtie bowtie2
~~~

### Install the gmap aligner

Install `gmap` using `conda`.

~~~ bash
conda install gmap
~~~

### Install the bwa aligner

Install `bwa` using `conda`.

~~~ bash
conda install bwa
~~~



# Installation of the ChimeraScan gene fusion prediction tool

## Environment setup

Set variable for index directory.

~~~ bash
CHIMERASCAN_INDEX=$TUTORIAL_HOME/refdata/chimerascan/indices/
~~~



## Installation

### Install ChimeraScan

Install in ChimeraScan using `conda`.

~~~ bash
conda install chimerascan
~~~

### Install the required reference data files

Install the reference data in a subdirectory of the tutorial ref data.

~~~ bash
mkdir -p $TUTORIAL_HOME/refdata/chimerascan/
cd $TUTORIAL_HOME/refdata/chimerascan/
~~~

Download the gene models from chimerascan's google code site as specified in the instructions.

~~~ bash
wget https://chimerascan.googlecode.com/files/hg19.ucsc_genes.txt.gz
gunzip hg19.ucsc_genes.txt.gz
~~~

Build the chimerascan indices using the `chimerascan_index.py` command.

~~~ bash
mkdir -p $CHIMERASCAN_INDEX
chimerascan_index.py \
    $UCSC_GENOME_FILENAME hg19.ucsc_genes.txt \
    $CHIMERASCAN_INDEX
~~~



# Installation of the deFuse gene fusion prediction tool

## Environment setup

Set variable for the config filename and the two scripts.

~~~ bash
DEFUSE_CONFIG=$TUTORIAL_HOME/cbw_tutorial/config/defuse_chr1.txt
DEFUSE_REF_DATA=$TUTORIAL_HOME/refdata/defuse/
~~~


## Installation

### Install deFuse

Install in ChimeraScan using `conda`.

~~~ bash
conda install defuse
~~~

### Install the required reference data files

The reference data files can be downloaded and index automatically using
the `defuse_create_ref.pl` script.

~~~ bash
defuse_create_ref.pl -c $DEFUSE_CONFIG -d $DEFUSE_REF_DATA
~~~



# Installation of the STAR RNA-Seq aligner

## Environment

Specify the directory in which the reference genome data will be stored.

~~~ bash
STAR_GENOME_INDEX=$TUTORIAL_HOME/refdata/star/
STAR_FUSION_GENOME_INDEX=$TUTORIAL_HOME/refdata/star/GRCh37_gencode_v19_CTAT_lib/
~~~

 

## Installation

Install STAR using conda.  Also install perl and the perl package Set::IntervalTree, required by STAR-Fusion.

~~~
conda install perl-threaded
conda install perl-set-intervaltree
conda install star
~~~

## Create genome 

For STAR-Fusion, we require an additional reference dataset.

~~~
cd $STAR_GENOME_INDEX

wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib.tar.gz

tar -xvf GRCh37_gencode_v19_CTAT_lib.tar.gz

cd GRCh37_gencode_v19_CTAT_lib/
~~~

Optionally subset the reference for the tutorial chromosomes.

~~~
samtools faidx ref_genome.fa chr$TUTORIAL_CHROMOSOME \
    > ref_genome.chr$TUTORIAL_CHROMOSOME.fa
mv ref_genome.chr$TUTORIAL_CHROMOSOME.fa ref_genome.fa
rm ref_genome.fa.fai

grep -P "^chr$TUTORIAL_CHROMOSOME\t" ref_annot.gtf \
    > ref_annot.chr$TUTORIAL_CHROMOSOME.gtf
mv ref_annot.chr$TUTORIAL_CHROMOSOME.gtf ref_annot.gtf
~~~

Prepare the genome and annotations for star fusion.

~~~
prep_genome_lib.pl \
    --genome_fa ref_genome.fa \
    --gtf ref_annot.gtf \
    --blast_pairs blast_pairs.outfmt6.gz
~~~




# Installation of the tophat-fusion gene fusion prediction tool

## Environment setup

Location of tophat-fusion specific gene models, ensembl but with chr prefix.

~~~ bash
TOPHAT_GTF_FILENAME=$TUTORIAL_HOME/refdata/tophatfusion/Homo_sapiens.GRCh37.75.gtf
~~~



## Installation

Install `tophat` and `tophat-fusion` using conda:

~~~ bash
conda install tophat
~~~

## Reference Genome Preparation

Install the reference gene annotations provided by tophat-fusion in the tophat-fusion specific
reference data directory.

~~~ bash
mkdir -p $TUTORIAL_HOME/refdata/tophatfusion/
cd $TUTORIAL_HOME/refdata/tophatfusion/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz
gunzip ensGene.txt.gz
~~~

We require the gene models in GTF format and we will use the GTF provided by ensembl.  However, we first need to add the 'chr' prefix to the chromosome names so we match the ucsc genome.  Use the sed command to create a modified version of the ensembl gene models with the chr prefix for each chromosome.  This can be done by just adding `chr` to the beginning of each line in the GTF
file.

~~~ bash
sed 's/^\([^#]\)/chr'$TUTORIAL_CHROMOSOME'/' $ENSEMBL_GTF_FILENAME | sed 's/^chrMT/chrM/' > $TOPHAT_GTF_FILENAME
~~~

Tophat requires additional blast databases, download and install these.

~~~ bash
mkdir -p $TUTORIAL_HOME/refdata/tophatfusion/blast
cd $TUTORIAL_HOME/refdata/tophatfusion/blast
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic.*.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz
gunzip *.gz
~~~




# Installation of the Trinity RNA-Seq assembler



## Installation

Install `trinity` using `conda`

~~~ bash
conda install trinity
~~~


