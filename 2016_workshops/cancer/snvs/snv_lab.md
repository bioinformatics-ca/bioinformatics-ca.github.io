# Lab Module 6 - Somatic Mutations

## Setup

First login into the server.

Now enter the ~/workspace directory

```
cd ~/workspace
```

Create a directory for this module and enter this directory:

```
mkdir Module6
cd Module6
```

Now let's create a link to some helper scripts we will need for this module:

```
ln -s /home/ubuntu/CourseData/CG_data/Module6/scripts
```

## Environment

In this section, we will set some environment variables to help facilitate the execution of commands. These variables will set to the location of some important files we need for the commands in this module. One important thing to remember is that:

> These variables that you set will only persist in this current session you are in. If you log out and log back into the server, you will have to set these variables again.

Set the directory of where MutationSeq is installed:

```
MUTATIONSEQ_DIR=/usr/local/museq/
```

Set the directory of where SnpEff is installed:

```
SNPEFF_DIR=/usr/local/snpEff
```

Summary of all the above environment commands (for copy and pasting convenience):

```
MUTATIONSEQ_DIR=/usr/local/museq
SNPEFF_DIR=/usr/local/snpEff
```

## Retrieving the Data

We will be using the same exome data on the  HCC1395 breast cancer cell line that was used in Module 5 (Copy Number Alterations).  The tumour and normal bam files have been placed on the server. Create a link to the location of the bam files to allow us to access them quickly.

```
ln -s /home/ubuntu/CourseData/CG_data/HCC1395
```

Additionally, we will need the reference genome for the mutation calling. This has also been placed on the server. We shall create a link to this file for easy access:

```
ln -s /home/ubuntu/CourseData/CG_data/ref_data
```

We will restrict our analysis to a 1 Mb region (7Mb and 8Mb) within chromosome 17. These bam files along with their indices have been generated for you already (please see the data preparation page for details on how to generate these smaller bams):

```
ls HCC1395/exome/HCC1395_exome*.17*
```

You should see:

* HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam
* HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam.bai
* HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam
* HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam.bai

If you are unfamiliar with the sam/bam format, take some time to look at the alignments and metadata.

The header information contains information about the reference genome and read groups (per read information about sample/flow cell/lane etc).

```
samtools view -H HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam | less -S
```

The main contents of the file contain read alignments, tab separated, one line per alignment.

```
samtools view HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam | less -S
```

Samtools will also calculate statistics about the reads and alignments.  Unfortunately this information is not cached, and thus this command will take considerable time on a regular sized bam.

```
samtools flagstat HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam
```

## Predicting SNVs

### Strelka

We will first call mutations using Strelka. Create a local copy of the Strelka config file.  Strelka provides aligner specific config files for bwa, eland, and isaac.  Each file contains default configuration parameters that work well with the aligner.  The bam files we are working with were created using bwa, so we select that config file and make a local copy to make changes.

```
mkdir config
cp /usr/local/etc/strelka_config_bwa_default.ini config/strelka_config_bwa.ini
```

Since we will be using exome data for this, we need to change the `isSkipDepthFilters` parameter in the strelka_config_bwa.ini file. Let's create a new config file for exome analysis:

```
cp config/strelka_config_bwa.ini config/strelka_config_bwa_exome.ini
```

Now let's edit the `config/strelka_config_bwa_exome.ini` and change the `isSkipDepthFilters = 0` to `isSkipDepthFilters = 1`.  We will use the text editor nano for this:

```
nano config/strelka_config_bwa_exome.ini
```

Please let us know if you have any issues editing the file. The reason why we do this is described on the [Strelka FAQ page](https://sites.google.com/site/strelkasomaticvariantcaller/home/faq):

> The depth filter is designed to filter out all variants which are called above a multiple of the mean chromosome depth, the default configuration is set to filter variants with a depth greater than 3x the chromosomal mean. If you are using exome/targeted sequencing data, the depth filter should be turned off...
> 
> However in whole exome sequencing data, the mean chromosomal depth will be extremely small, resulting in nearly all variants being (improperly) filtered out.

If you were doing this for whole genome sequencing data, then you should leave this parameter set to 0 as the depth of coverage won't be as high. 

A Strelka analysis is performed in 2 steps.  In the first step we provide Strelka with all the information it requires to run the analysis, including the tumour and normal bam filenames, the config, and the reference genome.  Strelka will create an output directory with the setup required to run the analysis.

```
configureStrelkaWorkflow.pl \
    --tumor HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam \
    --normal HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam \
    --ref ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --config config/strelka_config_bwa_exome.ini \
    --output-dir results/strelka/
```

The output directory will contain a _makefile_ that can be used with the tool _make_.  One benefit of the makefile style workflow is that it can be easily parallelized using _qmake_ on a grid engine cluster.  

To run the Strelka analysis, use make and specify the directory constructed by `configureStrelkaWorkflow.pl` with make's `'-C'` option.

```
make -C results/strelka/ -j 2
```

The `-j 2` parameter specifies that we want to use 2 cores to run Strelka. Change this number to increase or decrease the parallelization of the job. The more cores the faster the job will be, but the higher the load on the machine that is running Strelka. 

> If you have access to a grid engine cluster, you can replace the command `make` with `qmake` to launch Strelka on the cluster.

Strelka has the benefit of calling SNVs and small indels.  Additionally, Strelka calculates variant quality and filters the data in two tiers.  The filenames starting with `passed` contain high quality candidates, and filenames starting with `all` contain high quality and marginal quality candidates.

The Strelka results are in VCF format, with additional fields explained on the [strelka website](https://sites.google.com/site/strelkasomaticvariantcaller/home/somatic-variant-output).

```
less -S results/strelka/results/passed.somatic.snvs.vcf
less -S results/strelka/results/passed.somatic.indels.vcf
```

### MutationSeq

Now we will try another mutation caller called MutationSeq. MutationSeq uses supervised learning (random forest) to classify each mutation as true or artifact. The training set consists of over 1000 known true and positive mutations validated by deep SNV sequencing. The trained model is provided with the MutationSeq package, but must be provided on the command line using the `model:` argument.

Running MutationSeq is a one step process. 

```
mkdir -p results/mutationseq
python $MUTATIONSEQ_DIR/classify.py \
    tumour:HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam \
    normal:HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam \
    reference:ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    model:$MUTATIONSEQ_DIR/model_v4.1.1.npz \
    -i 17:7000000-8000000 \
    -c $MUTATIONSEQ_DIR/metadata.config -q 1 -o results/mutationseq/HCC1395.museq.vcf \
    -l results/mutationseq/HCC1395.log --all &
```

Some parameters to take note of:

- `-i` to specify the region
- `-q 1` to remove ambiguously mapping reads
- `--all` to print variants that are classified as low probability.
    
This job will run in the background. You can check its progress by going:

```
less -S results/mutationseq/HCC1395.log
```

When MutationSeq has finished, the results will be provided in VCF format:

```
less -S results/mutationseq/HCC1395.museq.vcf
```

Filtering on these results are left to the end-user. A threshold of 0.85 on PR (probability of being a true somatic mutation) will be used in this lab.

## Converting the VCF format into a tabular format

The VCF format is sometimes not useful for visualization and data exploration purposes which often requires the data to be in tabular format. We can convert from VCF format to tabular format using the extractField() function from SnpSift/SnpEff. Since each mutation caller has a different set of output values in the VCF file, the command needs be adjusted for the mutation caller. 

For example, to convert the Strelka VCF file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "."  results/strelka/results/passed.somatic.snvs.vcf CHROM POS ID REF ALT QUAL FILTER QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC GEN[0].DP GEN[1].DP GEN[0].FDP GEN[1].FDP GEN[0].SDP GEN[1].SDP GEN[0].SUBDP GEN[1].SUBDP GEN[0].AU GEN[1].AU GEN[0].CU GEN[1].CU GEN[0].GU GEN[1].GU GEN[0].TU GEN[1].TU > results/strelka/results/passed.somatic.snvs.txt 
```

To convert the MutationSeq VCF file into a tabular format:

```
java -jar $SNPEFF_DIR/SnpSift.jar extractFields -e "." results/mutationseq/HCC1395.museq.vcf CHROM POS ID REF ALT QUAL FILTER PR TR TA NR NA TC NI ND > results/mutationseq/HCC1395.museq.txt
```

The -e parameter specifies how to represent empty fields. In this case, the "." character is placed for any empty fields. This facilitates loading and completeness of data. For more details on the extractField() function see the [SnpSift documentation](http://snpeff.sourceforge.net/SnpSift.html#Extract).

## Data Exploration

### Integrative Genomics Viewer (IGV) 

A common step after prediction of SNVs is to visualize these mutations in IGV. Let's load these bam into IGV. Open IGV, then:

1. Change the genome to hg19 (if it isn't already)
2. File -> Load from URL ...
    * http://cbwxx.dyndns.info//Module6/HCC1395/exome/HCC1395_exome_tumour.17.7MB-8MB.bam
    * http://cbwxx.dyndns.info//Module6/HCC1395/exome/HCC1395_exome_normal.17.7MB-8MB.bam

Where the xx is your student number. Once the tumour and normal bam have been loaded into IGV, we can investigate a few predicted positions in IGV:

* 17:7491818
* 17:7578406
* 17:7482929

Manually inspecting these predicted SNVs in IGV is a good way to verify the predictions and also identify potential false positives:

> When possible, you should always inspect SNVs

### Exploration in R

While IGV is good for visualizing individual mutations, looking at more global characteristics would require loading the data into an analysis language like R:

We will use exome-wide SNV predictions for Strelka and MutationSeq for these analyses. These processed tabular text files along with the `analyzeSNVResults.Rmd` RMarkdown file that contains the R code for the analysis can downloaded as a package. 

Open the `analyzeSNVResults.Rmd` in RStudio now. 



