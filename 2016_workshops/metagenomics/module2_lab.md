---
layout: post2
permalink: /analysis_of_metagenomic_data_module2_lab_2016/
title: Analysis of Metagenomic Data 2016 Student Page
header1: Analysis of Metagenomic Data 2016
header2: Module 2 Lab
image: CBW_Metagenome_icon.jpg
---

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Module 2: Marker Gene Lab


Background
==========

In this lab, we will go over the major steps of 16S analysis using QIIME scripts and some additional custom scripts so we can become familiar with how to process and analyze 16S data. Analyze 16S data with QIIME entails running scripts each consist of multiple steps (commands). Some QIIME scripts even "call" (i.e. execute) additional QIIME scripts to form a workflow. As a result, running a single script could accomplish many tasks and generate multiple result files. This approach hides a lot of complexity and to really understand QIIME analyses, one has to study the documents available for the scripts (see http://qiime.org/scripts/index.html). This approach also make troubleshooting more difficult as one may have to understand the entire workflow to figure out which step is causing the error.

A example of a QIIME script is "split_libraries_fastq.py" and you can find the corresponding document at http://qiime.org/scripts/split_libraries_fastq.html. From the document, you can tell that the script performs demultiplexing of Fastq sequence data where barcodes and sequences are contained in two separate fastq files. To run the script, you would type "split_libraries_fastq.py -o slout/ -i forward_reads.fastq.gz -b barcodes.fastq.gz -m map.tsv" into your terminal window. Note that the script name is followed by a list of "agruments" (i.e. -o -i -b and -m).  These "agruments" (also called "parameters") change the default behavior of the script and provide additional information needed for the script to run.  For example, -i specifies the input sequence file to be processed and -b specifies the corresponding barcode file for the input sequences.  Information about these parameters is also provided in the script document page. Note that some arguments are optional and have default values (e.g. by default --store_demultiplexed_fastq is set to false so only a single file containing reads from all samples will be outputed). Lastly, the document has a section on the expected outputs of the script. 

For the lab, we will mainly use the dataset from the Mothur SOP (http://www.mothur.org/wiki/MiSeq_SOP) since it’s an interesting dataset and it has been pared down for demonstration purpose. Moreover, you could use the same dataset to run through the Mothur tutorial (bonus! at end of the QIIME tutorial) and learn to use both Mothur and QIIME. In this case, de-multiplexing has been done already so you have one pair of files (for paired-end reads) per sample. For Illumina MiSeq sequencing, de-multiplexing (or binning) may have been done for you already by the sequencing centre but if you need to do it on your own, both Mothur and QIIME have commands that you can use. If your data set has not been de-multiplexed (i.e. all your samples are in the same FASTQ file and you have the associated barcode file, then you can follow the QIIME Illumina tutorial (http://nbviewer.jupyter.org/github/biocore/qiime/blob/1.9.1/examples/ipynb/illumina_overview_tutorial.ipynb) to carry out the initial de-multiplexing and quality filtering of your data.

In this semi self-guided tutorial, you can copy and paste the commands in the grey boxes into your terminal to execute the commands. We will provide some explanation as you go along.  The learning objective here is to get you to be familiar with running interactive command line tools like QIIME.  In the Integrated Assignment, you will run all the commands from a single shell script - effectively running your entire analysis pipeline in one go and focus on interpreting the result files.


Dataset Intro (always good to know a bit about the data you are working with):
==============

The Schloss lab is interested in understanding the effect of normal variation in the gut microbiome on host health. Fresh feces from mice were collected for 365 days post weaning. In the first 150 days no intervention was done. In this demo, we will look at the data collected in the first 10 days except day 4 (early time points) vs. days 140-150 (late time points) for a single female mouse. In addition to mice fecal samples, we included a mock community composed of genomic DNA from 21 bacterial strains (available from HMP studies).

To connect to your Amazon Cloud Instance, type the following command in your terminal window (replace XX with your ID)

```
ssh -i CBW.pem ubuntu@cbwXX.dyndns.info
```


QIIME Workflow
==============

We will cover QIIME first rather than Mothur as it seemed to have gained more popularity over Mothur in the last few years.

In `~/CourseData/metagenomics/markergenes/qiime`, we have

**38 input sequences (FASTQ) files:**

-   F3Dxxx\_S209\_L001\_R1/2\_001.fastq.gz for the 19 mouse samples (forward and reverse reads)

**Reference Datasets:**

-   97_otus.fasta for GreenGenes reference database sequences (representative sequences from 97% identity cluster)
-   97_otu_taxonomy.txt for GreenGenes reference database taxonomy
-   core_set_aligned.fasta.imputed GreenGenes alignment template

**Metadata:** In QIIME, metadata is kept in a tab delimited file (.tsv) that's commonly called the mapping file

-   qiime_demo_metadata.tsv


First, we will setup our analysis directory.

```
mkdir -p ~/workspace/lab2_qiime
cd ~/workspace/lab2_qiime
```

Then we will link the files we need to perform the analysis. Linking files instead of copying the files to your analysis directory saves disk space. First make three directories to organize your input sequence files, your reference data, and the scripts that you'll be using.

```
mkdir sequence_files
mkdir reference_data
mkdir scripts
```

Now link the files from the ~/CourseData directory which is read-only.

```
ln -s ~/CourseData/metagenomics/markergenes/qiime/*.fastq.gz sequence_files/
ln -s ~/CourseData/metagenomics/markergenes/qiime/97* reference_data/
ln -s ~/CourseData/metagenomics/markergenes/qiime/core_set_aligned.fasta.imputed reference_data/
```

Then we will copy over the env file (metadata file) that describes our samples and a couple of scripts we need for QIIME analysis.

```
cp ~/CourseData/metagenomics/markergenes/qiime/qiime_demo_metadata.tsv ./
wget -O scripts/mesas-pcoa https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa
wget -O scripts/mesas-uc2clust https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-uc2clust
```

Then we make the scripts executable.

```
chmod u+x scripts/*
```

Lastly, we will setup our environmental variables so the scripts know where to find certain programs

```
export PYTHONPATH=~/local/lib/python2.7/site-packages
export RDP_JAR_PATH=/usr/local/rdp_classifier/classifier.jar
```


Pre-processing
--------------

### Paired-end Assembly

We will start with assemble our paired-end reads to merge the two reads into a single (possibly longer) read if the ends overlap. We will use a program called PEAR (Paired_End reAd mergeR) to do that. The document for PEAR can be found at http://sco.h-its.org/exelixis/web/software/pear/doc.html#cl-usage. While QIIME offers a script to assemble already de-multiplexed files (see http://qiime.org/tutorials/processing_illumina_data.html for more information), we will apply our own solution here to tell the assembly program (PEAR) how to match up the pairs of files and assemble the reads.


```
for i in "1" "5" "9" "13" "17" "21"
do
    echo $i
    find sequence_files/ -name "*.fastq.gz" -printf '%f\n' | sed 's/_L001.*//' | sort | uniq | sed -n $i,$((i+3))p | while read line; do ( pear -f sequence_files/${line}_L001_R1_001.fastq.gz -r sequence_files/${line}_L001_R2_001.fastq.gz -o ${line} & ); done > /dev/null
    sleep 40
done
```

While the above little snippet of code looks complicated, it basically tells the computer to interate through the sequence_files directory and run PEAR assembler on matching pairs of sequence files (XXX_L001_R1_001.fastq.gz and XXX_L001_R2_001.fastq.gz). Since there are 20 matching pairs, You could run the PEAR program yourself 20 times by changing the XXX to the correct file name (e.g. pear -f sequence_files/F3D0_S188_L001_R1_001.fastq.gz -r sequence_files/F3D0_S188_L001_R2_001.fastq.gz -o F3D0_S188). You can look at the files created by PEAR by typing

```
ls *.fastq
```

Notice that for each sample, there is one assembled read file, one discarded read file and two unassembled read files (forward and reverse).

As QIIME expects one single FASTA file as input sequence file, next we will combine all the assembled reads into a single file.  Again, this is done by looping through the directory containing the assembled.fastq files.  We also use a bit of "awk" magic to convert FASTQ files into a FASTA formatted file called "seq.fasta". Awk is a utility that allows one to process text files.  More information about the Awk language can be found at http://www.grymoire.com/Unix/Awk.html (come back to this later when you have time!)

```
for filename in $( ls *.assembled.fastq )
do
    awk 'BEGIN{ORS=""; i=0;}{split(FILENAME, x, "."); prefix=x[1]; sub("@","",prefix); print ">" prefix "_" i "\n"; i+=1; 
         getline; print; print "\n"; getline; getline;}' ${filename} >> seq.fasta
done
```

### Reduce the Number of Redundant Sequences

Next we will perform a set of procedures to reduce the number of redundant sequences and cluster our sequences (essentially forming OTUs). First we will collapse dataset to remove redundancy. We will do this using the default program for clustering in QIIME, namely usearch (note there is a new recommendation for clustering in QIIME. We will meet the new tools in the Integrated Assignment so you can compare the different clustering approaches).

```
usearch -derep_fulllength seq.fasta -fastaout derep.fa -sizeout
```

The result showed that from a total of 152132 sequences, there are 27,289 unique sequences.

00:00  83Mb  100.0% Reading seq.fasta
00:00  93Mb 152132 seqs, 27289 uniques, 22777 singletons (83.5%)
00:00  93Mb Min size 1, median 1, max 11708, avg 5.57
00:00  82Mb  100.0% Writing derep.fa


Then we will sort the sequences by abundance and cluster the sequences based on minimal 97% identity. In this step, chimera detection is also performed


```
usearch -sortbysize derep.fa -fastaout sorted.fa -minsize 2
usearch -cluster_otus sorted.fa -otus otus.fa -otu_radius_pct 3 -sizeout -uparseout results.txt
```


Here's the output. It shows that there are 219 OTUs based on the filtering and cutoffs we specified.
00:01  47Mb  100.0% 219 OTUs, 348 chimeras (7.7%)

Then to make the output QIIME compatible, we will need to rename the sequence IDs (to be numerical). Again, we rely on awk to do the text manipulation.

```
awk 'BEGIN{count=0;}{if ($$0~/>/){print ">" count; count+=1;} else {print}}' otus.fa>rep_set.fasta
```

This step, we will map all sequences (even singletons) back onto the OTUs. The file generated consists of OTU IDs and memberships. This file is similar in format to the .names. However, it is meant to store OTU membership.

```
usearch -usearch_global seq.fasta -db rep_set.fasta -strand both -id 0.97 -uc map.uc -threads 4
```


Lastly we will convert USEARCH output to QIIME OTU list format using a custom script (which you will have access to after the workshop)


```
scripts/mesas-uc2clust -t 4 map.uc seq_otus.txt
```


We will create a new directory to store the clusters we generated:


```
mkdir cluster
mv results.txt otus.fa map.uc rep_set.fasta seq_otus.txt sorted.fa derep.fa cluster
```


3: Taxonomy Assignment
----------------------

We will use Mothur Taxonomy (kmer based) Classifier for taxonomy assignment. Note that this is done by calling a QIIME script from command line directly. Like before, we need to provide a taxonomy file and a reference sequence file.

```
assign_taxonomy.py -m mothur -i cluster/rep_set.fasta -o taxonomic_classification -t reference_data/97_otu_taxonomy.txt -r reference_data/97_otus.fasta -c 0.6
```


The results of taxonomy assignment is in ~/workspace/lab2\_qiime/taxonomic\_classification/. As you can see the abundant OTUs are Bacteroidales. We use the sort command to sort the OTUs from most abundant to least abundant.


```
sort -k1n ./taxonomic_classification/rep_set_tax_assignments.txt | less
#less command allows you to scroll up and down the results
#type q to quick the scrolling window
```
The Results looks like:
```
0  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
1  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
2  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
3  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
4  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__ovatus 0.730
5  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
6  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__;s__   0.990
7  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
8  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
9  k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__S24-7;g__;s__   1.000
```

5: Sequence Alignment
---------------------

To align sequences to a template in order to do phylogenetic tree generation, we use align\_seqs.py. Unlike Mothur, the QIIME workflow does this after sequences are clustered as OTUs. The method (-m) used is pynast which is actually a very similar algorithm to what we used in Mothur. Pynast is a template based alignment algorithm so it requires a pre-aligned database (specified with -t option) to do the alignment.


```
align_seqs.py -m pynast -i cluster/rep_set.fasta -o alignment -t reference_data/core_set_aligned.fasta.imputed
filter_alignment.py -i alignment/rep_set_aligned.fasta -o alignment -s
```


6:Phylogenetic Analysis
-----------------------

We will first make a new directory to store the tree. Then we will run the make\_phylogeny.py script to create the tree from the aligned sequences using the fasttree algorithm. The fasttree algorithm as the name implied is really fast and can build tree from thousands of long sequences in a few minutes. It is much faster than the clearcut algorithm that mothur uses.


```
mkdir phylogeny
make_phylogeny.py -i alignment/rep_set_aligned_pfiltered.fasta -t fasttree -o phylogeny/rep_set.tre -l phylogeny/log.txt
```


4: Build OTU Table
------------------

We will first make a new directory to store the OTU table. Then we will run the script to make OTU table containing counts and taxonomic classifications


```
mkdir otu_table
make_otu_table.py -i cluster/seq_otus.txt -o otu_table/otu_table.biom -t taxonomic_classification/rep_set_tax_assignments.txt
```


This creates a .biom file in the otu\_table directory. The .biom file is a binary file so it is not readable directly. So we will convert the BIOM format file to tab-separated format (human readable). More information about the biom convert script can be found at <http://biom-format.org/documentation/biom_conversion.html>.


```
biom convert --to-tsv -i otu_table/otu_table.biom --table-type='OTU table' -o otu_table/otu_table.tab --header-key=taxonomy --output-metadata-id=Consensus\ Lineage
```


In order to carry out rarefaction of our data, we use awk on the converted OTU table to determine the lowest sequence depth. The result is passed as a parameter to QIIME's rarefaction script


```
single_rarefaction.py -i otu_table/otu_table.biom -o otu_table/otu_table_rarefied.biom -d `awk 'BEGIN{FS="\t"} NR == 1 { } NR == 2 { max = NF-1; } NR > 2 { for (i = 2; i <= max; i++) { c[i] += $i; } } END { smallest = c[2]; for (i = 3; i <= max; i++) { if (c[i] < smallest) { smallest = c[i]; }} print smallest; }' otu_table/otu_table.tab`
```
```
biom convert --to-tsv -i otu_table/otu_table_rarefied.biom --table-type='OTU table' -o otu_table/otu_table_rarefied.tab --header-key=taxonomy --output-metadata-id=Consensus\ Lineage
```

7:Downstream Analysis
---------------------

We can quickly plot the taxonomic profile of the samples by using the following command:

```
summarize_taxa_through_plots.py -f -s -i otu_table/otu_table.biom -o taxaplot -m qiime_demo_metadata.tsv
```


We can look at the output plots through a webpage (html file): <http://cbwXX.dyndns.info/lab2_qiime/taxaplot/taxa_summary_plots/bar_charts.html>. Remember to replace the XX with your machine number.

We will then calculate beta-diversity using weighted unifrac:


```
beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m weighted_unifrac -t phylogeny/rep_set.tre
```


Next we will generate the PCoA plot using these beta-diversity matrices


```
principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa/
make_2d_plots.py -i qiime_pcoa/pcoa -m qiime_demo_metadata.tsv -o qiime_pcoa/pcoa_plots/
chmod -R o+x qiime_pcoa/pcoa_plots
```


Lastly, we will run an R script for PCoA


```
mkdir mesas_pcoa
scripts/mesas-pcoa -i otu_table/otu_table_rarefied.tab -m qiime_demo_metadata.tsv -d bray -o mesas_pcoa/pcoa-bray.pdf
```


You can view the results at <http://cbwXX.dyndns.info/lab2_qiime/mesas_pcoa/pcoa-bray.pdf>. Again, remember to change the XX to your machine number.




Mothur Workflow (Optional)
===============

In `~/CourseData/markergenes/mothur`, we have

**40 input sequences (FASTQ) files:**

-   F3Dxxx\_S209\_L001\_R1/2\_001.fastq for the 19 mouse samples (forward and reverse reads)
-   Mock\_S280\_L001\_R1/2\_001.fastq for the mock community sample

**Reference Datasets:**

-   Silva (Silva): consists of both taxonomic files (.tax files) and representative reference sequence files (.fasta files)
-   RDP (trainset9): consists of Mothur formatted RDP 16S reference sequences and taxonomic assignments

**Metadata:** In Mothur, metadata are in two column formats (ID and Grouping)

-   stability.files (list of paired forward and reverse reads)
-   mouse.dpw.metadata (list of samples and day post weaning)
-   mouse.time.design (list of samples and early or late time points)

Adapted from [Mothur MiSeq SOP](http://www.mothur.org/wiki/MiSeq_SOP). For full tutorial please visit the Mothur website.

First, we will create our workspace for today's demo.


```
mkdir -p ~/workspace/lab2_mothur
cd ~/workspace/lab2_mothur
```


After you are in the ~/workspace/lab2_mothur directory, we will copy over the necessary files.  Linking is probably a better option


```
cp ~/CourseData/markergenes/mothur/* ~/workspace/lab2_mothur/
```


Then we will start mothur by typing "mothur"


```
mothur
```


**Since mothur is going to occupy your current terminal, we will open another terminal to be able to view the files generated by mothur commands.**

1: Pre-processing
-----------------

### Paired-end Assembly

Because we have overlapping paired-end reads, one way to improve the overall sequence quality is by assembling the paired-end reads. The command to do that in Mothur is


```
mothur > make.contigs(file=stability.files, processors=8)
```


Mothur uses its own simple algorithm to merge the paired reads. We will see later in QIIME workflow a different tool (pear) for merging paired-end reads. The command output a bunch of messages to the screen, but what's worth paying attention is the list of output files generated by mothur


```
Output File Names:  
stability.trim.contigs.fasta
stability.contigs.qual
stability.contigs.report
stability.scrap.contigs.fasta
stability.scrap.contigs.qual
stability.contigs.groups
```


Mothur will always give you a list of files it generated - a handy feature to have. Note that all the reads are now put into a single fasta file called "stability.trim.contigs.fasta". The information linking reads to the samples are in "stability.contigs.groups". We can take a quick look at the assembly results by:


```
mothur > summary.seqs(fasta=stability.trim.contigs.fasta)
```


This command which we will see over and over again summarizes the content of a fasta sequence file. We can see a total of 152,360 sequences successfully assembled. We can also notice that most reads are within the expected size (~250bps) but some reads are way too long. These reads are poorly assembled. We will further improve the input sequence quality by issuing the following command:


```
mothur > screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=275)
```


This command remove sequences that are too long or contain ambiguous bases (e.g. N's). Again, we will run the summary.seqs command to check the output.


```
mothur > summary.seqs(fasta=stability.trim.contigs.good.fasta)
```


You will notice that we now have 128,872 that passes the filter.

### Reduce the Number of Redundant Sequences

In your samples, there will be sequences that are identical. We can issue the following commands to combine duplicated sequences into one and keep track of the number of copies of each of the duplicated sequences.


```
mothur > unique.seqs(fasta=stability.trim.contigs.good.fasta)
mothur > count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
```


This will generate "stability.trim.contigs.good.unique.fasta" which contains only one copy of each of the duplicated reads. We can see there are only 16426 unique sequences. The "stability.trim.contigs.good.names" file keeps track of the sequences that are duplicated. The count.seqs command generate a count\_table file "stability.trim.contigs.good.count\_table" that just keep track of the number of copies of each of the duplicated sequences. In general, the .names file gives you the membership of each cluster of sequences.

### Alignment Sequences to Reference Database

Next, we will align the sequences to the silva reference alignment.


```
mothur > align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.bacteria.fasta)
```


We'll see the alignment summary with:


```
mothur > summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table)
```


You can see most of the sequences align to region between 13862 and 23444. We will remove poorly aligned sequences and remove gap-only columns by issuing the following commands:


```
mothur > screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=13862, end=23444, maxhomop=8)
mothur > filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.)
```


The alignment and filter results are:


```
Length of filtered alignment: 375
Number of columns removed: 49625
Length of the original alignment: 50000
Number of sequences used to construct filter: 16300
```


After alignment and trimming, some of the reads may now be duplicated, so we will run the unique.seqs again.


```
mothur > unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table)
```


We now have aligned sequences in "stability.trim.contigs.good.unique.good.filter.unique.fasta" and the duplicated sequences are tracked in "stability.trim.contigs.good.unique.good.filter.count\_table".

We will use pre.cluster command to further remove rare sequences that are very similar to abundant sequences (within 2nt differences). These rare sequences are likely to be due to sequencing errors.


```
mothur > pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2)
```


### Remove Chimeric and non-16S sequences

Next we will use Mothur to remove chimeric sequences.


```
mothur > chimera.uchime(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
mothur > remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos)
```


Using the classify.seqs command (which is mothur's version of RDP Classifier), we can remove sequences that do not look like a 16S sequence from bacteria. Note this step requires the RDP database (trainset9\_032012.pds.fasta and the taxonomy file) as reference.


```
mothur > classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)
mothur > remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
```


Note that the classify.seqs step is akin to closed-reference OTU picking. The sequences are matched the a reference sequence in the database and unmatched sequences can be removed using remove.lineage or remove.seqs commands.

2: OTU Picking
--------------

Now we are done with pre-processing, the next stage of our analysis is to make OTUs.

Mothur takes a de-novo OTU picking approach by first calculate the pair-wise distances of each of the sequences then cluster them into OTUs using hierarchical clustering. Notice that by specifying a cutoff of 0.20, only pairs of sequences that are 80% identical will be recorded. Since we will not be interested in clustering together dissimilar sequneces, this cutoff will help us save disk and memory spaces. Also notice that we need the aligned sequences to calculate the distance matrix.


```
mothur > dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, cutoff=0.20)
mothur > cluster(column=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table)
```


The cluster step generate "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique\_list.list" file which contain the membership of each OTUs.

Next we want to know how many sequences are in each OTU from each group and we can do this using the make.shared command. Here we tell mothur that we're really only interested in the 0.03 cutoff level:


```
mothur > make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, label=0.03)
```


This creates a .shared file ([Mothur:SharedFile](http://www.mothur.org/wiki/Shared_file)) which contains the OTU count (similar content as the OTU table in QIIME).

If we open the .shared file, you can see there are 402 OTUs.

3: Taxonomy Assignment
----------------------

In Mothur, taxonomy assignment can be done using the classify.otu command which takes the OTU list and count\_table (to keep track of the read counts in each OTU) and searches against the taxonomy results we got from classify.seqs to assign OTUs to taxons.


```
mothur > classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, label=0.03)
```


You can see the results of taxonomy assignment in "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique\_list.0.03.cons.taxonomy". Go to the other terminal window and issue the following command:


```
less stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy
```



```
OTU     Size    Taxonomy
Otu0001 12328   Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100); 
Otu0002 8915    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0003 7850    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0004 7479    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0005 7478    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);Barnesiella(100);
Otu0006 6650    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Porphyromonadaceae"(100);unclassified(100);
Otu0007 6341    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);Bacteroidaceae(100);Bacteroides(100);
Otu0008 5374    Bacteria(100);"Bacteroidetes"(100);"Bacteroidia"(100);"Bacteroidales"(100);"Rikenellaceae"(100);Alistipes(100);
```


From the results, you can see Porphyromonadaceae are the most common organisms found in our samples.

4: Build OTU Table
------------------

In Mothur, the OTU Tables are represented by the .shared files. We can convert .shared files to QIIME compatible OTU tables in BIOM format with the following command: (see [Make.biom](http://www.mothur.org/wiki/Make.biom))


```
mothur > make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy)
```


5: Sequence Alignment
---------------------

In Mothur, the sequence alignment step is done earlier as part of the pre-processing. This allowed Mothur to perform de-novo OTU picking on the sequence data.

6: Phylogenetic Analysis
------------------------

If you are interested in using methods that depend on a phylogenetic tree such as calculating phylogenetic diversity or the unifrac commands, you'll need to generate a tree. In mothur, clearcut is used to generate a phylogenetic tree of the OTUs.


```
mothur > dist.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, output=lt, processors=8)
mothur > clearcut(phylip=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.dist)
```


This generates a .tre file (stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre). You can download this file by going to <http://cbwXX.dyndns.info/lab2_mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre> (replacing the XX with your machine number)

You can open the .tre file using [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) program.

7: Downstream Analysis
----------------------

Now we have our taxonomy file and our OTU table, we will quickly rename these files to make the downstream analysis a bit more legible


```
mothur > system(mv stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared stability.an.shared)
mothur > system(mv stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy stability.an.cons.taxonomy)
```


We will then Rarefy our samples so each sample has the same number of reads:


```
mothur > count.groups(shared=stability.an.shared)
mothur > sub.sample(shared=stability.an.shared, size=2440)
```


### Alpha Diversity

We can look at the diversity within each sample by generating the rarefaction curve:


```
mothur > rarefaction.single(shared=stability.an.shared, calc=sobs, freq=100)
```


This will generate files ending in \*.rarefaction, which again can be plotted in your favorite graphing software package (let's try with Excel). Alas, rarefaction is not a measure of richness, but a measure of diversity. If you consider two communities with the same richness, but different evenness then after sampling a large number of individuals their rarefaction curves will asymptote to the same value. Since they have different evennesses the shapes of the curves will differ.

We can generate a summary of the alpha-diversity by:


```
mothur > summary.single(shared=stability.an.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=2440)
```


### Beta Diversity

Beta diversity measures the differences in communities memberships or community structures across samples: We will first use the dist.shared function to calculate the distances between pairs of samples. ThetaYC is a measure of community structure and Jclass is a measure of community membership (i.e. presence or absence of OTUs)


```
mothur > dist.shared(shared=stability.an.shared, calc=thetayc-jclass, subsample=2240)
```


Then we will turn the distance calculation into a bifurcating tree.


```
mothur > tree.shared(phylip=stability.an.thetayc.0.03.lt.ave.dist)
```


We can visualize the tree in FigTree again. From the tree, we can see that the early and late samples cluster separately.

We can also take phylogenetic information into account when measuring beta-diversity. This is done by using UniFrac


```
mothur > unifrac.unweighted(tree=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, distance=lt, processors=2, random=F, subsample=2240)
mothur > unifrac.weighted(tree=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, distance=lt, processors=2, random=F, subsample=2240)
```


These commands will distance matrices (stability.phylip.tre1.unweighted.ave.dist and stability.phylip.tre1.weighted.ave.dist) that can be analyzed using all of the beta diversity approaches described above for the OTU-based analyses. For example, we could use the following command to generate a sample bifurcating tree.


```
 mothur > tree.shared(phylip=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.phylip.tre1.unweighted.phylip.dist)
```


There are other downstream analyses such as picking out indicator species and classifying community types. I will leave these to your own exploration.

To quit Mothur, simply type

```
mothur > quit
```


Summary
=======

We ran the same dataset through fairly similar Mothur and QIIME pipelines. While the exact number of OTUs are not the same (probably due to the inclusion and exclusion of singletons and the average linkage vs. minimal linkage clustering approaches), the resulting patterns from diversity plots and PCoA are consistent. The lab is to give you a taste of how these tools work, you can explore the results more on your own and during the integrated assignment. Feel free to explore and ask questions!
