# Lab Module 5 - Data Preparation

## Introduction

As mentioned in the lab, the data and some of the pre-processing steps have been done for you already to save time. In this page, we provide the instructions on how to retrieve such data and how the pre-processing was done. Use the information on this page for reference purposes in case you wish to reproduce the lab in your own time. 

## Environment

In this section, we will set some environment variables to help facilitate the execution of commands. These variables will set to the location of some important files we need for the commands in this module. One important thing to remember is that:

> These variables that you set will only persist in this current session you are in. If you log out and log back into the server, you will have to set these variables again.

Set the directory of where SnpEff is installed:

```
SNPEFF_DIR=/path/to/SnpEff
```

Set the directory of where HMMCopy is installed:

```
HMMCOPY_DIR=/path/to/hmmCopy
```

Set the directory of where Samtools is installed:

```
SAMTOOLS_DIR=/path/to/samtools
```

The pre-processing helper scripts can be downloaded from the wiki. Set the directory of where these scripts are stored:

```
CNA_SCRIPTS_DIR=/path/to/scripts
```

## Analysis Of CNAs using Arrays

In this section, we will cover the retrieval and pre-processing of some of the data used in the lab section 
"Analysis of CNAs using Arrays"

### Downloading the HCC1395 Affymetrix 6.0 Microarray

The CEL files containing the probe intensity data for cell line HCC1395 can be downloaded from the NCBI Geo website as follows:

```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM888nnn/GSM888107/suppl/GSM888107.CEL.gz
gunzip GSM888107.CEL.gz
```

## Analysis Of CNAs using Sequencing Data

In this section, we cover the retrieval and pre-processing of some of the data used in the lab section 
"Analysis Of CNAs using Sequencing Data"

### Downloading the HCC1395 Whole Exome Data

The exome data along with the other sequencing data can be downloaded from `https://github.com/genome/gms/wiki/HCC1395-WGS-Exome-RNA-Seq-Data`

### Calculate Tumour and Normal Read Depth using HMMCopy

Copy number prediction is based on read depth.  We can generate this information using HMMCopy to calculate read depth for the tumour and normal bam files and store this information in `.wig` format.  We use the -c parameter to specify all autosome chromosomes:

```
mkdir -p hmmCopy/wig 
$HMMCOPY_DIR/bin/readCounter -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 HCC1395/exome/HCC1395_exome_tumour.bam \ 
	> hmmCopy/wig/HC1395_exome_tumour.wig
$HMMCOPY_DIR/bin/readCounter -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 HCC1395/exome/HCC1395_exome_normal.bam \
	> hmmCopy/wig/HC1395_exome_normal.wig 
```

### Retrieve Allele Counts of Heterozygous Positions from the Tumour and Normal

The next step we need to do is retrieve the allele count data from heterozygous positions in both the tumour and normal. The first we will do is get a list of all heterozygous positions in the normal:

```
mkdir -p titan/bcftools/vcf 
$SAMTOOLS_DIR/samtools mpileup -u -I -f ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa HCC1395/exome/HCC1395_exome_normal.bam | \
	$SAMTOOLS_DIR/bcftools/bcftools view -vcg - | \
	java -jar $SNPEFF_DIR/SnpSift.jar filter "isHet(GEN[0]) & (QUAL >= 20)" > \
	titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf
```

Here we used samtools and bcftools (1.18) to predict the genotype at each position followed by SnpSift to select the heterozygous positions. We impose a `QUAL >= 20` filter that only retrieves positions with a quality of > 20 to enhance quality. You can use other programs (e.g. GATK) to retrieve these heterozygous positions. Another way to enhance quality is to intersect your list of positions with reported positions in dbSNP. We do not do this in this case, but SnpSift provides functionality to intersect your positions with a list of positions.

Once we have the normal heterozygous positions, we need to find the corresponding information at these positions in the tumour. The first thing we do is generate a `bed` file of these positions:

```
mkdir -p titan/bcftools/bed
grep -v "^#" titan/bcftools/vcf/HCC1395_exome_normal.var.het.vcf | awk '{OFS="\t"; print $1,$2-1,$2,$4"/"$5,"+"}' > titan/bcftools/bed/HCC1395_exome_normal.var.het.bed
```

This `bed` file provides a way to pass as input the positions we want to retrieve data from in the tumour:

```
$SAMTOOLS_DIR/samtools view -uF 1024 HCC1395/exome/HCC1395_exome_tumour.bam | \
	$SAMTOOLS_DIR/samtools mpileup -u -I -f ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -l titan/bcftools/bed/HCC1395_exome_normal.var.het.bed - | \
	$SAMTOOLS_DIR/bcftools/bcftools view -cg - > titan/bcftools/vcf/HCC1395_exome_tumour.bcftools.vcf
```

Now we extract the relevant information out of the tumour VCF file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields titan/bcftools/vcf/HCC1395_exome_tumour.bcftools.vcf CHROM POS REF ALT DP4 > titan/bcftools/tables/HCC1395_exome_normal.var.het.table.tmp
```

And finally, we convert this tabular format into the exact input format that TITAN needs. We use a perl script for this:

```
perl $CNA_SCRIPTS_DIR/getTitanAlleleReadCounts.pl titan/bcftools/tables/HCC1395_exome_normal.var.het.table.tmp
 > titan/bcftools/tables/HCC1395_exome_normal.var.het.table.txt
```

### Generating the Genome Reference Mappability File

We use HMMCopy to generate the genome reference mappability file. This a two step process where we first need to generate a `BigWig` file followed by converting it into a `Wig` file

```
$HMMCOPY_DIR/util/mappability/generateMap.pl refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw
```

This may fail if you have not built a bowtie index of your reference yet. To do this, make sure first that bowtie is in your $PATH variable. Once this is done, you can go:

```
$HMMCOPY_DIR/util/mappability/generateMap.pl -b refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa -o refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw
```

This will build the bowtie index. Now you can re-run the first command. In the end, you should get a BigWig file which you then need to convert into a wig file. You can do this by going:

```
$HMMCOPY_DIR/bin/mapCounter -w 1000 refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.bw > refdata/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig
```

### Generating the Genome Reference GC Content File

We use HMMCopy to generate the genome reference GC content file. This is one step command:

```
$HMMCOPY_DIR/bin/gcCounter Homo_sapiens.GRCh37.75.dna.primary_assembly.fa > Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig
```

## CNA Data Analysis Pre-processing

In this section, we will cover how some of the output data from OncoSNP and TITAN was pre-processed for the data analysis portion of the lab.