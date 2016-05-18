---
layout: post3
title: Bioinformatics for Cancer Genomics 2016 Gene Fusions Tutorial
header1: Bioinformatics for Cancer Genomics 2016
header2: Gene Fusions Tutorial
image: fusions.png
---

# General Setup

## Create tutorial directory structure on the instance

Create a directory in the temporary workspace and change to that directory.

~~~ bash
mkdir /mnt/workspace/Module4/
cd /mnt/workspace/Module4/
~~~

Gene fusion prediction tools have been installed on the instance, create
links to the installed binaries, libraries, and reference genome data.

~~~ bash
ln -s /media/cbwdata/CG_data/Module4/tutorial_data/bin
ln -s /media/cbwdata/CG_data/Module4/tutorial_data/lib/
ln -s /media/cbwdata/CG_data/Module4/tutorial_data/packages
ln -s /media/cbwdata/CG_data/Module4/tutorial_data/refdata
ln -s /media/cbwdata/CG_data/Module4/tutorial_data/share
~~~

Set the `$TUTORIAL_HOME` environment variable so subsequent commands know the 
locatioin of installed binaries, libraries, and reference genome data.  Also,
subsequent analyses will fall into subdirectories of `$TUTORIAL_HOME`.

~~~ bash
TUTORIAL_HOME=/mnt/workspace/Module4
~~~


## Environment setup

For this tutorial, we now assume the environment variable `$TUTORIAL_HOME`
has been set to an existing directory to which the user has write access.

Add tutorial binaries to the path.

~~~ bash
export PATH=$TUTORIAL_HOME/bin:$PATH
~~~

Create a variable pointing to the directory of the picard tools install.

~~~ bash
PICARD_DIR=$TUTORIAL_HOME/packages/picard-tools-1.130/
~~~

Create a variable pointing to the directory of the igvtools install.

~~~ bash
IGVTOOLS_DIR=$TUTORIAL_HOME/packages/IGVTools/
~~~



## Install additional libraries

### R ada for deFuse

Install the additional R package 'ada' locally.

~~~ bash
mkdir /mnt/workspace/Module4/rlib/

Rscript -e "install.packages('ada', \
	repos='http://cran.us.r-project.org', \
	lib='/mnt/workspace/Module4/rlib/')"
~~~

Add the local R library path to the `$R_LIBS_USER` environment variable.
Use `export` so that it is available to subprocesses, in particular the
`R` subprocess.

~~~ bash
export R_LIBS_USER=/mnt/workspace/Module4/rlib/
~~~

### `Set::IntervalTree` for STAR-Fusion

Install the `Set::IntervalTree` perl module required by star-fusion.

~~~ bash
mkdir $TUTORIAL_HOME/perlpackages
cd $TUTORIAL_HOME/perlpackages

git clone https://github.com/amcpherson/Set-IntervalTree.git
cd Set-IntervalTree

mkdir $TUTORIAL_HOME/perllib
perl Makefile.PL PREFIX=$TUTORIAL_HOME/perllib LIB=$TUTORIAL_HOME/perllib
make
make install
~~~

Add the local perl library path to the `$PERL5LIB` environment variable.
Use `export` so that it is available to perl subprocesses.

~~~ bash
export PERL5LIB=$PERL5LIB:$TUTORIAL_HOME/perllib
~~~

### Python `jinja2` for chimerascan

First install the `jinja2` python package locally.

~~~ bash
pip install jinja2 -t /mnt/workspace/Module4/pythonlib
~~~

Add the install directory to the python path environment variable so
that python can find the package.

~~~ bash
export PYTHONPATH=$PYTHONPATH:$TUTORIAL_HOME/pythonlib
~~~





## Environment setup

Create variables for the reference genome and gene models.

~~~ bash
TUTORIAL_CHROMOSOME=1
ENSEMBL_GTF_FILENAME=$TUTORIAL_HOME/refdata/Homo_sapiens.GRCh37.75.gtf
GENCODE_GTF_FILENAME=$TUTORIAL_HOME/refdata/gencode.v19.annotation.gtf
ENSEMBL_GENOME_FILENAME=$TUTORIAL_HOME/refdata/Homo_sapiens.GRCh37.75.dna.chromosome.$TUTORIAL_CHROMOSOME.fa
UCSC_GENOME_FILENAME=$TUTORIAL_HOME/refdata/chr$TUTORIAL_CHROMOSOME.fa
ENSEMBL_GENOME_PREFIX=$TUTORIAL_HOME/refdata/Homo_sapiens.GRCh37.75.dna.chromosome.$TUTORIAL_CHROMOSOME
UCSC_GENOME_PREFIX=$TUTORIAL_HOME/refdata/chr$TUTORIAL_CHROMOSOME
~~~

For gmap, set a directory for the gmap indices, and the id of the reference.

~~~ bash
GMAP_REF_ID=chr$TUTORIAL_CHROMOSOME
GMAP_INDEX_DIR=$TUTORIAL_HOME/refdata/gmap/
~~~



# HCC1395 sample data

## Environment setup

We will create environment variables pointing to the data on the instance, and
set a sample id for the data.

We will be working with a publically available HCC1395 breast cancer cell line
sequenced for teaching and benchmarking purposes.  The dataset is available 
from washington university servers, and can be accessed via the github page for the
[Genome Modeling System](https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data).

Special thanks to Malachi Griffith for providing access to the data.

~~~ bash
SAMPLE_ID=HCC1395
SAMPLE_FASTQ_1=/media/cbwdata/CG_data/HCC1395/rnaseq/SAMPLE_R1.fastq
SAMPLE_FASTQ_2=/media/cbwdata/CG_data/HCC1395/rnaseq/SAMPLE_R2.fastq
~~~



# Running times for each tool

chimerascan: 12 minutes

defuse: 18 minutes

star: 11 minutes

tophatfusion: 17 minutes

trinity: 12 minutes


# Run the ChimeraScan fusion prediction tool

## Environment setup

Add local library to python path

~~~ bash
export PYTHONPATH=$PYTHONPATH:$TUTORIAL_HOME/lib/python2.7/site-packages
~~~

Set variable for index directory.

~~~ bash
CHIMERASCAN_INDEX=$TUTORIAL_HOME/refdata/chimerascan/indices/
~~~



## Execution

Run chimerascan by providing the index directory, the pair of fastq files, and
the sample specific output directory

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/chimerascan/$SAMPLE_ID/
python $TUTORIAL_HOME/bin/chimerascan_run.py -p 4 \
    $CHIMERASCAN_INDEX \
    $SAMPLE_FASTQ_1 $SAMPLE_FASTQ_2 \
    $TUTORIAL_HOME/analysis/chimerascan/$SAMPLE_ID/
~~~

After chimerascan has completed we can run an additional script to create an
html output to browse the results.

~~~ bash
python $TUTORIAL_HOME/bin/chimerascan_html_table.py --read-throughs \
	-o $TUTORIAL_HOME/analysis/chimerascan/$SAMPLE_ID/chimeras.html \
	$TUTORIAL_HOME/analysis/chimerascan/$SAMPLE_ID/chimeras.bedpe
~~~



# Run the deFuse gene fusion prediction tool

## Environment setup

Set variable for the config filename and the two scripts.

~~~ bash
DEFUSE_CONFIG=$TUTORIAL_HOME/refdata/defuse/config.txt
CREATEREF_SCRIPT=$TUTORIAL_HOME/packages/defuse/scripts/create_reference_dataset.pl
DEFUSE_SCRIPT=$TUTORIAL_HOME/packages/defuse/scripts/defuse.pl
DEFUSE_SCRIPT_DIR=$TUTORIAL_HOME/packages/defuse/scripts/
~~~


## Execution

Running deFuse involves invoking a single script using perl.  The output is to
a directory which will contain a number of temporary and results files.  Create
an output directory and run the defuse script specifying the paired end reads
and the configuration filename.

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/defuse/$SAMPLE_ID/

perl $DEFUSE_SCRIPT \
	-c $DEFUSE_CONFIG \
	-1 $SAMPLE_FASTQ_1 \
	-2 $SAMPLE_FASTQ_2 \
	-o $TUTORIAL_HOME/analysis/defuse/$SAMPLE_ID/
~~~


# Run the STAR RNA-Seq aligner on a sample

## Environment

Specify the directory in which the reference genome data will be stored.

~~~ bash
STAR_GENOME_INDEX=$TUTORIAL_HOME/refdata/star/
~~~

Create a variable for the location of the star-fusion perl script for easy running of star-fusion.

~~~ bash
STAR_FUSION_SCRIPT=$TUTORIAL_HOME/packages/STAR-Fusion/STAR-Fusion
~~~



## Execution

### Run the star aligner

Running STAR on paired end read fastqs can be performed in a single step.  Set the prefix
argument `--outFileNamePrefix` to specify the location of the output.

Assume we have end 1 and end 2 files as `$SAMPLE_FASTQ_1` and `$SAMPLE_FASTQ_2` environment
variables pointing to paired fastq files for sample `$SAMPLE_ID`.

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/
STAR --runThreadN 4 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMstrandField intronMotif \
	--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
	--outReadsUnmapped None \
	--chimSegmentMin 15 \
	--chimJunctionOverhangMin 15 \
	--chimOutType WithinBAM \
    --alignMatesGapMax 200000 \
    --alignIntronMax 200000 \
	--genomeDir $STAR_GENOME_INDEX \
	--readFilesIn $SAMPLE_FASTQ_1 $SAMPLE_FASTQ_2 \
	--outFileNamePrefix $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/
~~~

> Note: In order to obtain alignments of chimeric reads potentially supporting fusions, we
> have added the `--chimSegmentMin 20` option to obtain chimerica reads anchored by at least
> 20nt on either side of the fusion boundary, and `--chimOutTypeWithinBAM` to report such
> alignments in the sam/bam output.

> Note: We are running with additional command line options for fusion discovery as described
> in the manual for STAR-Fusion.

The file `$TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Aligned.sortedByCoord.out.bam` should contain
the resulting alignments in bam format, sorted by position.

Index the bam file using `samtools index`.

~~~ bash
samtools index $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Aligned.sortedByCoord.out.bam
~~~

Generate a precalculated coverage file for use with IGV.  This will 
allow viewing of read depth across the genome at all scales, and will
speed up IGV viewing of the bam file.

~~~ bash
$IGVTOOLS_DIR/igvtools count $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Aligned.sortedByCoord.out.bam \
	$TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Aligned.sortedByCoord.out.bam.tdf hg19
~~~

You can view the bam file in IGV at:

~~~
http://cbw#.dyndns.info/Module4/analysis/star/HCC1395/Aligned.sortedByCoord.out.bam
~~~

where you need to replace # with your student ID.

### Run STAR fusion caller

The STAR fusion caller parses the chimeric alignments produced by the STAR aligner
and predicts fusions from these alignments.

~~~ bash
perl $STAR_FUSION_SCRIPT \
	--chimeric_out_sam $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.sam \
	--chimeric_junction $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.junction \
	--ref_GTF $GENCODE_GTF_FILENAME \
	--out_prefix $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/starfusion
~~~

The fusion reads are output to the sam file `Chimeric.out.sam`.  To view these in
IGV, first import into bam format, sort, and then index.

~~~ bash
samtools view -bt $UCSC_GENOME_FILENAME \
	$TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.sam \
	| samtools sort - $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out
samtools index $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.bam
$IGVTOOLS_DIR/igvtools count $TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.bam \
	$TUTORIAL_HOME/analysis/star/$SAMPLE_ID/Chimeric.out.bam.tdf hg19
~~~

You can view the bam file of fusion reads in IGV at:

~~~
http://cbw#.dyndns.info/Module4/analysis/star/HCC1395/Chimeric.out.bam
~~~

where you need to replace # with your student ID.  Fusion reads can be found
at this location `chr1:176,901,821-176,904,807`.




# Run the tophat-fusion gene fusion prediction tool

## Environment setup

Set the package dependent on the system, linux or osx.

~~~ bash
if [[ `uname` == 'Linux' ]]; then
    TOPHAT_PACKAGE=tophat-2.0.14.Linux_x86_64
elif [[ `uname` == 'Darwin' ]]; then
    TOPHAT_PACKAGE=tophat-2.0.14.OSX_x86_64
fi
~~~

Location of tophat-fusion specific gene models, ensembl but with chr prefix.

~~~ bash
TOPHAT_GTF_FILENAME=$TUTORIAL_HOME/refdata/tophatfusion/Homo_sapiens.GRCh37.75.gtf
~~~



## Execution

Running Tophat-fusion is a two step process.  First we run tophat2 to do spliced alignments, adding
the `--fusion-search` option to ensure we are looking for fusion splice alignments.  Remember to use
the prefix of the reference genome, rather than the fasta filename.

Run the tophat2 step specifying a sample specific output directory.

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/tophatfusion/$SAMPLE_ID/

tophat2 -p 4 \
    --mate-inner-dist 100 \
    --mate-std-dev 30 \
    --fusion-search \
    -G $TOPHAT_GTF_FILENAME \
    -o $TUTORIAL_HOME/analysis/tophatfusion/$SAMPLE_ID/ \
    $UCSC_GENOME_PREFIX \
    $SAMPLE_FASTQ_1 $SAMPLE_FASTQ_2
~~~

The second step requires us to soft link some of the reference data into a directory with
a specific structure.  We then run the tophat-fusion post-processing step.

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/tophatfusion/$SAMPLE_ID/tophat_$SAMPLE_ID
cd $TUTORIAL_HOME/analysis/tophatfusion/$SAMPLE_ID/tophat_$SAMPLE_ID

ln -s $TUTORIAL_HOME/refdata/tophatfusion/blast
ln -s $TUTORIAL_HOME/refdata/tophatfusion/refGene.txt
ln -s $TUTORIAL_HOME/refdata/tophatfusion/ensGene.txt

tophat-fusion-post \
    --num-fusion-reads 1 \
    --num-fusion-pairs 2 \
    --num-fusion-both 5 \
    $UCSC_GENOME_PREFIX
~~~



# Run the Trinity RNA-Seq assembler / GMAP contig mapper

## Environment

The trinity package contains some perl modules, add the directory containing these modules
to the perl library directory, so perl knows where to find them.

~~~ bash
export PERL5LIB=$TUTORIAL_HOME/packages/trinityrnaseq_r20140413p1/PerlLib/:$PERL5LIB
~~~



## Execution

### Assemble using Trinity

Create a directory for the trinity analysis.  Run trinity with fastq (fq) as the 
sequence type.  Specify memory and threads, and provide the fastq paths and output
paths.

~~~ bash
mkdir -p $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/

cp $SAMPLE_FASTQ_1 $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/sample_reads_1.fq
cp $SAMPLE_FASTQ_2 $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/sample_reads_2.fq

Trinity \
	--seqType fq --JM 12G --CPU 4 \
	--left $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/sample_reads_1.fq \
	--right $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/sample_reads_2.fq \
	--output $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/
~~~

### Align using GMap

Use gmap to map the resulting contigs to the reference genome and produce a GFF3
file.

~~~ bash
cd $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/

gmap \
	-f gff3_gene \
	-D $GMAP_INDEX_DIR \
	-d $GMAP_REF_ID \
	Trinity.fasta \
	> Trinity.gff3
~~~

### Postprocess

Post-process the gff3 file to produce a list of fusions.  We will use a custom
script to simply pull out chimeric contigs.  Realistically, further processing
would be required to identify true fusions.

~~~ bash
cd $TUTORIAL_HOME/analysis/trinity/$SAMPLE_ID/

python $TUTORIAL_HOME/packages/cbw_tutorial/gene_fusions/scripts/gmap_extract_fusions.py \
	Trinity.gff3 Trinity.fasta Fusions.fasta
~~~

The following sequence is the ASTN1-FAM5B fusion.

~~~
>c107_g1_i2 len=396 path=[1:0-245 247:246-395]
GTTTAGAGGGCTTCGGCCGGGGATGGCCCCATGGACAGCCCTGCTGGCACTGGGCCTGCC
TGGCTGGGTGTTGGCTGTCTCAGCCACGGCGGCTGCTGTGGTCCCCGAGCAGCATGCCTC
CGTAGCTGGCCAGCATCCCCTGGACTGGCTGCTCACAGACCGGGGCCCCTTCCACCGCGC
TCAGGAGTATGCTGACTTCATGGAGCGGTACCGCCAGGGTTTCACCACCAGGTACAGGAT
TTATAGCCAGAGTCGCCGCTGGACCTGCTTGTTTGGGTACCAGATAGAACAGGACTCCTC
AAAGCCGTAGATAGCTTCCTGGATGTAATGGTTGCCGAACTGGTCCAACAGCGCCACAAA
ATCTGCACGAGATGTAGCCCCATCCAGCGAGTGGAG
~~~




