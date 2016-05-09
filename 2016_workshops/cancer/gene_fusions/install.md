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

Add tutorial binaries to the path.

```{.bash}
export PATH=$TUTORIAL_HOME/bin:$PATH
```

Create a variable pointing to the directory of the picard tools install.

```{.bash}
PICARD_DIR=$TUTORIAL_HOME/packages/picard-tools-1.130/
```

Create a variable pointing to the directory of the igvtools install.

```{.bash}
IGVTOOLS_DIR=$TUTORIAL_HOME/packages/IGVTools/
```

 

## Installation

### Tutorial directory structure

Create directories for packages, binaries, reference data.

```{.bash}
mkdir -p $TUTORIAL_HOME/packages
mkdir -p $TUTORIAL_HOME/bin
mkdir -p $TUTORIAL_HOME/refdata
```

### Install tutorial scripts

Clone the tutorial repo so we have a copy of the tutorial scripts in a known location.

```{.bash}
cd $TUTORIAL_HOME/packages
```

```{.bash}
git clone https://bitbucket.org/dranew/cbw_tutorial.git
```

### Install samtools

The `samtools` package is the most widely used software for manipulating high-throughput
sequence data stored in the 'bam' format.

Download the source to the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

```{.bash}
wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2 -O samtools-1.2.tar.bz2
tar -xjvf samtools-1.2.tar.bz2
```

Build samtools and install to the `$TUTORIAL_HOME/bin` directory.

```{.bash}
cd samtools-1.2
make
make prefix=$TUTORIAL_HOME install
```

### Install picard tools

Picard tools is a useful set of utilities for manipulating sequence data in bam/sam format.
To run picard tools you need:

* [Java Runtime Environment](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)

With the above installed, download the picard tools java archives.

```{.bash}
cd $TUTORIAL_HOME/packages
wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip -O picard-tools-1.130.zip
unzip picard-tools-1.130.zip
```

### Install igv tools

The igvtools package provides utilities for preprocessing bam files for quicker viewing in IGV.

```{.bash}
cd $TUTORIAL_HOME/packages
wget http://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.52.zip
unzip igvtools_2.3.52.zip
```

### Install bowtie and bowtie2

Install bowtie dependent on platform.

```{.bash}
cd $TUTORIAL_HOME/packages
if [[ `uname` == 'Linux' ]]; then
    BOWTIE_PACKAGE=bowtie-1.1.1-linux-x86_64
elif [[ `uname` == 'Darwin' ]]; then
    BOWTIE_PACKAGE=bowtie-1.1.1-macos-x86_64
fi
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.1/$BOWTIE_PACKAGE.zip
unzip $BOWTIE_PACKAGE.zip
cp bowtie-1.1.1/bowtie* $TUTORIAL_HOME/bin/
```

Install bowtie2 dependent on platform.

```{.bash}
cd $TUTORIAL_HOME/packages
if [[ `uname` == 'Linux' ]]; then
    BOWTIE2_PACKAGE=bowtie2-2.2.5-linux-x86_64
elif [[ `uname` == 'Darwin' ]]; then
    BOWTIE2_PACKAGE=bowtie2-2.2.5-macos-x86_64
fi
wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/$BOWTIE2_PACKAGE.zip
unzip $BOWTIE2_PACKAGE.zip
cp bowtie2-2.2.5/bowtie2* $TUTORIAL_HOME/bin/
```

### Install the gmap aligner

Install in the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

Download and unpack the source.

```{.bash}
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-12-29.tar.gz
tar -xzvf gmap-gsnap-2014-12-29.tar.gz
```

Build the source.

```{.bash}
cd gmap-2014-12-29
./configure --prefix=$TUTORIAL_HOME
make
make install
```

### Install the bwa aligner

Install in the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

Download and unpack the source.

```{.bash}
wget --no-check-certificate http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2
tar -xjvf bwa-0.7.12.tar.bz2
```

Build the source.

```{.bash}
cd bwa-0.7.12
make
cp bwa $TUTORIAL_HOME/bin/
```



# Installation of the ChimeraScan gene fusion prediction tool

## Environment setup

Add local library to python path

```{.bash}
export PYTHONPATH=$PYTHONPATH:$TUTORIAL_HOME/lib/python2.7/site-packages
```

Set variable for index directory.

```{.bash}
CHIMERASCAN_INDEX=$TUTORIAL_HOME/refdata/chimerascan/indices/
```



## Installation

### Install the source code

Install in the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

Download the chimerascan code.

```{.bash}
wget --no-check-certificate https://chimerascan.googlecode.com/files/chimerascan-0.4.5a.tar.gz
tar -xzvf chimerascan-0.4.5a.tar.gz
```

Enter the chimerascan directory to build source.

```{.bash}
cd chimerascan-0.4.5
```

Use sed for a minor patch for OSX compilation.

```{.bash}
sed 's/inline/static inline/' chimerascan/pysam/samtools/ksort.h
```

```{.bash}
python setup.py install --prefix $TUTORIAL_HOME
```

### Install the required reference data files

Install the reference data in a subdirectory of the tutorial ref data.

```{.bash}
mkdir -p $TUTORIAL_HOME/refdata/chimerascan/
cd $TUTORIAL_HOME/refdata/chimerascan/
```

Download the gene models from chimerascan's google code site as specified in the instructions.

```{.bash}
wget https://chimerascan.googlecode.com/files/hg19.ucsc_genes.txt.gz
gunzip hg19.ucsc_genes.txt.gz
```

Build the chimerascan indices using the `chimerascan_index.py` command.

```{.bash}
mkdir -p $CHIMERASCAN_INDEX
python $TUTORIAL_HOME/bin/chimerascan_index.py \
    $UCSC_GENOME_FILENAME hg19.ucsc_genes.txt \
    $CHIMERASCAN_INDEX
```



# Installation of the deFuse gene fusion prediction tool

## Environment setup

Set variable for the config filename and the two scripts.

```{.bash}
DEFUSE_CONFIG=$TUTORIAL_HOME/refdata/defuse/config.txt
CREATEREF_SCRIPT=$TUTORIAL_HOME/packages/defuse/scripts/create_reference_dataset.pl
DEFUSE_SCRIPT=$TUTORIAL_HOME/packages/defuse/scripts/defuse.pl
DEFUSE_SCRIPT_DIR=$TUTORIAL_HOME/packages/defuse/scripts/
```


## Installation

### Install the source code

Clone the defuse code into the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
git clone http://bitbucket.org/dranew/defuse.git
```

Change to the tools directory to build the tools from source.

```{.bash}
cd $TUTORIAL_HOME/packages/defuse/tools
make
```

Create a copy of the default config and modify some of the entries.  We will
do this using sed, but equivalently you could use your favourite text editor.

```{.bash}
mkdir $TUTORIAL_HOME/refdata/defuse/

cat $TUTORIAL_HOME/packages/defuse/scripts/config.txt \
    | sed 's#\[Where you unpacked the defuse code\]#'$TUTORIAL_HOME/packages/defuse'#' \
    | sed 's#\[Where you intend to store the reference dataset\]#'$TUTORIAL_HOME'/refdata/defuse#' \
    | sed 's#\[path of your \(.*\) binary.*\]#\1#' \
    > $DEFUSE_CONFIG
```

### Install the required reference data files

The reference data files can be downloaded and index automatically using
the `create_reference_dataset.pl` script.

```{.bash}
perl $CREATEREF_SCRIPT -c $DEFUSE_CONFIG
```



# Installation of the STAR RNA-Seq aligner

## Environment

Specify the directory in which the reference genome data will be stored.

```{.bash}
STAR_GENOME_INDEX=$TUTORIAL_HOME/refdata/star/
```

Create a variable for the location of the star-fusion perl script for easy running of star-fusion.

```{.bash}
STAR_FUSION_SCRIPT=$TUTORIAL_HOME/packages/STAR-Fusion/STAR-Fusion
```

 

## Installation

Clone the star repo in the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

```{.bash}
git clone --recursive https://github.com/alexdobin/STAR.git
```

Copy the star executable.

```{.bash}
if [[ `uname` == 'Linux' ]]; then
    cp STAR/bin/Linux_x86_64_static/STAR $TUTORIAL_HOME/bin
elif [[ `uname` == 'Darwin' ]]; then
    cp STAR/bin/MacOSX_x86_64/STAR $TUTORIAL_HOME/bin
fi
```

Clone the star-fusion repo and checkout version 0.1.1.

```{.bash}
cd $TUTORIAL_HOME/packages
```

```{.bash}
git clone https://github.com/STAR-Fusion/STAR-Fusion.git
cd STAR-Fusion
git checkout v0.1.1
```

Install the `Set::IntervalTree` perl module required by star-fusion.

```{.bash}
cd $TUTORIAL_HOME/packages
```

```{.bash}
mkdir $TUTORIAL_HOME/lib
git clone https://github.com/amcpherson/Set-IntervalTree.git
cd Set-IntervalTree
perl Makefile.PL PREFIX=$TUTORIAL_HOME/lib LIB=$TUTORIAL_HOME/lib
make
make install
```

## Create genome 

Index the genome.  Here we need to include the reference genome fasta and the gene models in GTF format.

```{.bash}
mkdir $STAR_GENOME_INDEX
STAR --runThreadN 1 \
    --runMode genomeGenerate \
    --genomeDir $STAR_GENOME_INDEX \
    --genomeFastaFiles $UCSC_GENOME_FILENAME \
    --sjdbGTFfile $GENCODE_GTF_FILENAME \
    --sjdbOverhang 100
```



# Installation of the tophat-fusion gene fusion prediction tool

## Environment setup

Set the package dependent on the system, linux or osx.

```{.bash}
if [[ `uname` == 'Linux' ]]; then
    TOPHAT_PACKAGE=tophat-2.0.14.Linux_x86_64
elif [[ `uname` == 'Darwin' ]]; then
    TOPHAT_PACKAGE=tophat-2.0.14.OSX_x86_64
fi
```

Location of tophat-fusion specific gene models, ensembl but with chr prefix.

```{.bash}
TOPHAT_GTF_FILENAME=$TUTORIAL_HOME/refdata/tophatfusion/Homo_sapiens.GRCh37.75.gtf
```



## Installation

Install in the packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

Download and extract the source code.

```{.bash}
wget --no-check-certificate wget --no-check-certificate http://ccb.jhu.edu/software/tophat/downloads/$TOPHAT_PACKAGE.tar.gz
tar -xzvf $TOPHAT_PACKAGE.tar.gz
```

Build the source code and install binaries.

```{.bash}
cd $TOPHAT_PACKAGE
ln -s `pwd`/tophat2 $TUTORIAL_HOME/bin/tophat2
ln -s `pwd`/tophat-fusion-post $TUTORIAL_HOME/bin/tophat-fusion-post
```

## Reference Genome Preparation

Install the reference gene annotations provided by tophat-fusion in the tophat-fusion specific
reference data directory.

```{.bash}
mkdir -p $TUTORIAL_HOME/refdata/tophatfusion/
cd $TUTORIAL_HOME/refdata/tophatfusion/
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
gunzip refGene.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ensGene.txt.gz
gunzip ensGene.txt.gz
```

We require the gene models in GTF format and we will use the GTF provided by ensembl.  However,
we first need to add the 'chr' prefix to the chromosome names so we match the ucsc genome.  Use
the sed command to create a modified version of the ensembl gene models with the chr prefix for
each chromosome.  This can be done by just adding `chr` to the beginning of each line in the GTF
file.

```{.bash}
sed 's/^\([^#]\)/chr\1/' $ENSEMBL_GTF_FILENAME | sed 's/^chrMT/chrM/' > $TOPHAT_GTF_FILENAME
```

Tophat requires additional blast databases, download and install these

```{.bash}
mkdir -p $TUTORIAL_HOME/refdata/tophatfusion/blast
cd $TUTORIAL_HOME/refdata/tophatfusion/blast
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/human_genomic.*.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz
gunzip *.gz
```




# Installation of the Trinity RNA-Seq assembler

## Environment

The trinity package contains some perl modules, add the directory containing these modules
to the perl library directory, so perl knows where to find them.

```{.bash}
export PERL5LIB=$TUTORIAL_HOME/packages/trinityrnaseq_r20140413p1/PerlLib/:$PERL5LIB
```



## Installation

Install into tutorial packages directory.

```{.bash}
cd $TUTORIAL_HOME/packages
```

Download and extract the trinity source.

```{.bash}
wget --no-check-certificate http://sourceforge.net/projects/trinityrnaseq/files/trinityrnaseq_r20140413p1.tar.gz
tar -xzvf trinityrnaseq_r20140413p1.tar.gz
```

Make the executable and add a symbolic link in the bin directory.

```{.bash}
cd trinityrnaseq_r20140413p1
make
ln -s `pwd`/Trinity $TUTORIAL_HOME/bin/Trinity
```


