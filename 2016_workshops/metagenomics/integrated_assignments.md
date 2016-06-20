---
layout: post2
permalink: /analysis_of_metagenomic_data_integrated_assignments_2016/
title: Analysis of Metagenomic Data 2016 Student Page
header1: Analysis of Metagenomic Data 2016
header2: Integrated Assignments
image: CBW_Metagenome_icon.jpg
---

<ul id="navmenu">
  <li><a id="back_to_top">Contents</a>
     <ul class="sub1">
     <li><a href="#background">Background</a></li>
      <li><a href="#day1">Day 1</a>
         <ul class="sub2">  
           <li><a href="#setup">Setup</a></li>
           <li><a href="#workflow">Workflow</a></li>
           <li><a href="#output">Output File Glossary</a></li>
           <li><a href="#questions">Questions</a></li>
        </ul>
      </li>
       <li><a href="#day_2">Day 2</a>
          <ul class="sub2">
             <li><a href="#who">Who is Out There</a></li>
             <li><a href="#what">What are They Doing</a></li>
             <li><a href="#questions2">Questions</a></li>
           </ul>
       </li>
    </ul>
  </li>
</ul>  

<br>

Background <a id="background"></a>
==========

This integrated lab assignment session will use a subset of samples collected and sequenced through the MicroB3 consortium (http://www.microb3.eu/osd) as part of the global Ocean Sampling Day event in 2014.

To quote from the OSD website:

   "The Ocean Sampling Day (OSD) is a simultaneous sampling campaign of the world’s oceans which took place (for the first time) on the summer solstice (June 21st) in the year 2014. These cumulative samples, related in time, space and environmental parameters, provide insights into fundamental rules describing microbial diversity and function and contribute to the blue economy through the identification of novel, ocean-derived biotechnologies. We see OSD data as a reference data set for generations of experiments to follow in the coming decade."

We will be analyzing samples collected from marine stations from various locations in the North West Atlantic shelf and the Polar Atlantic Arctic province. For the integrated assignment on day 1, we will be using sequences from the V4-V6 region of the 16S rRNA gene from 28 samples (19 from the North West Atlantic shelf and 9 from the Polar Atlantic Arctic). These samples have been subsampled from the originals so that the analysis runs faster. On day 2, we will be using shotgun metagenomic sequences from 25 samples (16 from the North West Atlantic shelf and 9 from the Polar Atlantic Arctic). These samples have been downloaded for you already from the European Bioinformatics Institute metagenomics portal (https://www.ebi.ac.uk/metagenomics/projects/ERP009703) and placed in the directories:

```
~/CourseData/integrated\_assignment\_day1/

~/CourseData/integrated\_assignment\_day2/
```

Objective 
---------

We have chosen a subset of the OSD samples from 2 regions which are expected to differ from each other in terms of temperature and probably also in terms of anthropogenic impact with the samples from the North-west Atlantic shelf being more impacted by human activity than the Arctic samples. The idea is to investigate if these differences are also reflected in a change in the microbial communities in these samples and if so, what are the microbial taxa or metabolic processes that drive these differences.

Integrated Assignment - Day \#1 <a id="day1"></a>
===============================

Set up <a id="setup"></a>
------

Once you have logged into the server, you will create the workspace for today's assignment.

```
mkdir -p ~/workspace/assignment1
cd ~/workspace/assignment1
```

All of the required commands are in a BASH script which we can copy to our workspaces. Each of the commands contained within this script are displayed and discussed below.

```
cp ~/CourseData/metagenomics/integrated_assignment_day1/Integrated_Lab_1.sh .
chmod u+x Integrated_Lab_1.sh
```

The first command fetches the script and the second command allows us to execute it.

This script takes about 40-60 minutes to run to completion. When running the script, piping to the `tee` command is helpful to store a copy of the output for later viewing.

```
./Integrated_Lab_1.sh 2>&1 | tee -a log.txt
```

Run the commands above then come back to this document to learn about what it's doing.

**You should not enter the commands below into the console. They are being run using the script.**

Let's go through the script command by command. First is a bunch of setup commands.

```
#Setup some parameters for the script
dataLocation="~/CourseData/metagenomics/integrated_assignment_day1/" # this is where your fastq files should be
workingDir="~/workspace/assignment1"
ncores=4

cd $workingDir

#First, we will link our sequence and reference data into our workspace
mkdir sequence_files
mkdir reference_data
mkdir scripts
ln -s $dataLocation/*.fastq.gz sequence_files/
ln -s $dataLocation/nwcs_arct_metadata.tsv .
ln -s $dataLocation/97* reference_data/
ln -s $dataLocation/*.py scripts
ln -s $dataLocation/mesas-* scripts # from [https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa](https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa)
```

These commands get the necessary data files into our workspace. Source FASTQ sequence files are linked to the sequence\_files directory. GreenGenes reference set and alignment template are linked to the reference\_data directory.

```
#Prep database for SortMeRNA using: cd reference_data; /usr/local/sortmerna-2.1/indexdb_rna --ref 97_otus.fasta
sortmernaDB=$workingDir"/reference_data/97_otus"

#Mark the scripts as executable
chmod u+x scripts/*
```

Here the `chmod` command allows us to execute the scripts.

Workflow <a id="workflow"></a>
--------

### Data QC

The first step when analysing data is to assess its quality and to check your assumptions. Are the reads of the expected length? Are the PHRED quality scores sufficiently high? We can run common quality checks using a tool called FastQC. This is a general quality checker for all types of DNA sequencing data, so not all of the output will be important for amplicon reads. Focus on the read lengths and quality scores.

```
fastqc -q -t $ncores sequence_files/*.fastq -o fastqc_out/raw/individuals
```

This command will run FastQC on all of our raw data files and store the results in the folder indicated. $ncores is a variable defined at the begining of the script indicating how many cores we have, here we're using that number to tell FastQC how many threads to use.

```
gunzip --to-stdout sequence_files/*.fastq.gz | fastqc -q -t $ncores stdin -o fastqc_out/raw/combined
```

If you have hundreds of data files, it can be overwhelming to check all the quality scores. If your data is \*perfect\*, you may not need to look at all the files individually. The command above runs FastQC on all the reads together and puts them in one output. If the results look \*perfect\*, then you're ready to go (best to always look at your data, however). If they don't look perfect, you'll need to go through the output from the previous command to find out which samples are problematic. A slow but crucial step!

### Paired-end Assembly

The next step in our pipeline is to assemble the Illumina paired-end reads with [PEAR](https://github.com/xflouris/PEAR).

```
ncores=4`
numberOfFilePairs=$(ls $workingDir/sequence_files/*_1.fastq.gz| wc -l )
let numberOfIterations=($numberOfFilePairs+$ncores-1)/$ncores
for j in $(seq 0 $numberOfIterations)
do
   let i=( $j * $ncores + 1 )
   echo $i
   find $workingDir/sequence_files -name "*.fastq.gz" -printf '%f\n' | sed 's/_.*//' | sort | uniq | sed -n $i,$((i+${ncores}-1))p | while read line; do ( pear -f sequence_file$
   sleep 60
done
```

Don't panic! This complex command is meant to showcase the power of the tools available in most BASH shell environments. This command is finding the forward and reverse sequence read files and launching 4 simultaneous PEAR processes every minute. By launching multiple PEAR processes at once, we save processing time. The complicated run-up with `find`, `sed`, `sort`, and `uniq` finds the files we want to assemble automatically so that we don't have to type all 54 file names in manually. A simpler approach would be to copy and paste the PEAR command, replacing the forward and reverse file arguments each time.

```
while [ $( ps -ef | grep pear | wc -l ) -gt 1 ]
do
 sleep 10
 echo "Waiting on PEAR to finish..."
done
```

This command detects if there are any PEAR processes running, and if there are, it will wait for them to finish.

```
grep "^Assembled reads" logPear.txt > logPearAssemblyStats.txt
```

This command gets the percentage of reads assembled from the PEAR log files. You should check through this result file manually to make sure a high percentage of the pairs assembled. If the percentages are low, or there is a wide range of percentages, or there is one sample with an unusually low percentage, you should follow up by looking at FastQC results of the pre-merged files to figure out why this is the case before using any results.

It's also important to keep track of what's happening to your data. In the previous step we merged read pairs and in the next step we're only going to use those reads that were successfully merged - effectively discarding reads that were left unmerged. If there's a problem with your data, you may end up discarding lots of data here and severely biasing your results (we've seen it before!).

```
for filename in $( ls *.assembled.fastq )
do
    awk 'BEGIN{ORS=""; i=0;}{split(FILENAME, x, "."); prefix=x[1]; sub("@","",prefix); print ">" prefix "_" i "\n"; i+=1; 
         getline; print; print "\n"; getline; getline;}' ${filename} >> seq.fasta
done
```

Here we leverage `awk` to convert the assembled FASTQ files into FASTA files where the sequence headers are in the format "&gt;SAMPLEID\_SEQUENCENUMBER". This is the format that QIIME expects for its input sequence files.

```
mkdir sequence_files
mv *.fastq sequence_files
```

Staying organized is important! It is easy to create hundreds of files and forget which is which.

### Generating OTUs with QIIME: Classification, Alignment, Trees, and OTU Tables

Next, we start the QIIME clustering pipeline. QIIME is a software package of python wrapper scripts which makes it easy to cluster sequences into OTUs, align OTU representative sequences to make a tree, and taxonomically classify OTUs.

We will use "open reference clustering", which first clusters our sequences against a database of 16S references sequences, then uses *de novo* clustering on those sequences which were not similar to the reference sequences. For more detail on this method and other options, see the [QIIME documentation](http://qiime.org/1.9.0/tutorials/otu_picking.html).

As QIIME is mostly a (very effective) package of wrapper scripts calling other software, there are choices to be made about what underlying software to use at every step of the pipeline. For example, the default clustering software is UCLUST, but in this case we will use a combination of SortMeRNA and SUMACLUST, as recommended in [this paper](http://msystems.asm.org/content/1/1/e00003-15) from the lab that develops QIIME.

You can see the options for the open reference clustering by looking at the [script documentation](http://qiime.org/scripts/pick_open_reference_otus.html). There is similar documentation for [other QIIME scripts](http://qiime.org/1.9.0/scripts/index.html) with many more details on each step of the pipeline.

Let's get back to reading through our script:

```
inputFasta=sequence_files/combined_fasta/combined_seqs.fna
```

This is the fasta file we generated in the previous section. We will use it as the input for OTU generation.

First, we set up parameters for how we want QIIME to run and save them in a file which we will pass to QIIME. For more on parameter files, see the [documentation](http://qiime.org/documentation/qiime_parameters_files.html).

```
echo "pick_otus:threads " $ncores >> clustering_params.txt
echo "pick_otus:sortmerna_coverage 0.8" >> clustering_params.txt
echo "assign_taxonomy:id_to_taxonomy_fp $workingDir/reference_data/97_otu_taxonomy.txt" >> clustering_params.txt
echo "pick_otus:sortmerna_db $sortmernaDB" >> clustering_params.txt
```

This last line is pointing to the SortMeRNA database, which was created from the 97\_otus.fasta reference file using this command: `sortmerna-2.1/indexdb\_rna --ref 97\_otus.fasta`. If you want to use a different reference database, you'll need to create a new database using this command.

```
pick_open_reference_otus.py -i $inputFasta -o clustering/ -p clustering_params.txt -m sortmerna_sumaclust -s 0.1 -v --min_otu_size 1
```

This is a QIIME wrapper script which will generate OTUs from our data and save the results in a directory called "clustering". Check out details on the script [here](http://qiime.org/1.9.0/scripts/pick_open_reference_otus.html), where the six major steps are outlines. This QIIME script calls other QIIME scripts for each of these steps, including some which are responsible for:

-   clustering sequences into OTUs (pick\_otus.py)
  -   
-   taxonomically classifying OTUs using the (assign\_taxonomy.py). The default setting is to use RDP against the [GreenGenes 13\_8 reference](ftp://greengenes.microbio.me/greengenes_release/unversioned/).
-   aligning the OTU representative sequences and filtering the alignment of columns that are all gaps (align\_seqs.py, filter\_alignment.py)
-   generating a tree of OTUs from the alignment, e.g. useful for UniFrac analysis (make\_phylogeny.py)

This script produces many results files, the most commonly used include:

-   otu\_table\_mc1\_w\_tax\_no\_pynast\_failures.biom : the final OTU results, including taxonomic assignments and per-sample abundances, stored in a biom file. This is the file you will use the most.
-   rep\_set.tre : a tree file of OTUs in newick format, indicating the similarity-based hierarchical relationship between OTUs
-   rep\_set.fna : a fasta file of the OTU representative sequences
-   final\_otu\_map\_mc1.txt : listing of which reads were clustered into which OTUs
-   pynast\_aligned\_seqs/rep\_set\_aligned\_pfiltered.fasta : alignment of OTU representative sequences
-   uclust\_assigned\_taxonomy/rep\_set\_tax\_assignments.txt : taxonomic assignment for each OTU

```
scripts/remove_low_confidence_otus.py -i clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o clustering/otu_table_high_conf.biom
```

Next, we use a custom script from [Microbiome Helper](https://github.com/mlangill/microbiome_helper/wiki) (led by Dr. Morgan Langille) that will discard low abundance OTUs. These may be due to cross-over from MiSeq lanes or sequencing errors OR they may be real biological low abundance taxonomic groups. Unless your study focuses on low abundance taxa, it's best to discard these extremely low abundance OTUs to facilitate downstream analysis and avoid technical artefacts. Make sure this is appropriate for your study before running this on your data, though!

```
biom summarize-table -i clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o clustering/otu_table_mc1_w_tax_no_pynast_failures_summary.txt
biom summarize-table -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_summary.txt
```

Results from clustering are stored in a "biom" file. Here we prepare summaries that describe the biom files before and after we discarded low abundance OTUs. When the script has finished, try running the `head` command with the two summaries and note the difference in the number of reads and the number of OTUs in each.

```
biom convert --to-tsv -i clustering/otu_table_high_conf.biom --table-type='OTU table' -o clustering/otu_table_high_conf.tsv --header-key=taxonomy --output-metadata-id=Consensus\ Lineage
```

Here we convert the BIOM format file to tab-separated format (human readable) using a command from the "biom" [software package](http://biom-format.org/) with a column for taxonomic info called "Consensus Lineage".

```
#We use awk on the converted OTU table to determine the lowest sequence depth
subsampleSize=$(awk 'BEGIN{FS="\t"} NR == 1 { } NR == 2 { max = NF-1; } NR > 2 { for (i = 2; i <= max; i++) { c[i] += $i; } } \
 END { smallest = c[2]; for (i = 3; i <= max; i++) { if (c[i] < smallest) { smallest = c[i]; }} print smallest; }' clustering/otu_table_high_conf.tsv)
echo $subsampleSize

#This is passed as a parameter to QIIME's rarefaction script
single_rarefaction.py -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_rarefied.biom -d $subsampleSize

```

The awk command uses the converted OTU table to determine the lowest sequence depth of all the input samples and we then pass that value as the -d parameter to QIIME's rarefaction script. This QIIME script will "rarefy" a biom file such that every sample in the output biom file (here: clustering/otu\_table\_high\_conf\_rarefied.biom) will be represented by the same number of reads.  This number could also be read manually from the summary files we generated.

### Downstream Analyses

We have our OTU tables! All our hard work has paid off, and we can now explore the data with visualizations and other analyses.

```
biomFile=clustering/otu_table_high_conf_rarefied.biom
biomTable=clustering/otu_table_high_conf_rarefied.tsv

summarize_taxa_through_plots.py -f -s -i $biomFile -o taxaplot -m nwcs_arct_metadata.tsv
```

Here we create bar plots of OTUs coloured by taxa.

```
beta_diversity_through_plots.py -m nwcs_arct_metadata.tsv -t clustering/rep_set.tre -i $biomFile -o qiime_pcoa_3D/
```

Next we create some 3D PCoA plots. You can view them by transferring them to your local computer or by going to [<http://cbwXX.dyndns.info>](http://cbwXX.dyndns.info) (where XX is your student number). These plots are useful to get a sense of the data but not so good for publications. The next few commands use QIIME scripts to give us 2D plots. (Alternatively, you can bring the biom file into R and use the phyloseq and/or vegan packages.)

```
beta_diversity.py -i $biomFile -o qiime_pcoa/distance/ -m weighted_unifrac -t clustering/rep_set.tre
beta_diversity.py -i $biomFile -o qiime_pcoa/distance/ -m bray_curtis
```

The first step is to create the distance matrices. Note that for UniFrac, a tree is required since it is a phylogenetic beta diversity measure.

```
#Run PCoA on these matrices
principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa/
```

Next, we run the PCoA algorithm on the data matrices (this is a QIIME script). This will rotate and scale the data in order to capture as much variance as possible in two (or three) dimensions.

```
make_2d_plots.py -i qiime_pcoa/pcoa -m nwcs_arct_metadata.tsv -o qiime_pcoa/pcoa_plots/
chmod -R o+x qiime_pcoa/pcoa_plots
```

QIIME will also plot this for us, and create a unique colouring for each metadata category in the nwcs\_arct\_metadata.tsv mapping file.

You can look at the results from these analyses using: [<http://cbwXX.dyndns.info>](http://cbwXX.dyndns.info) (where XX is your student number).


If you will be using R to analyse amplicon data, you will likely find the following packages useful: 

- [phyloseq](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)
- [vegan](http://cran.r-project.org/web/packages/vegan/index.html)
- [ape](http://cran.r-project.org/web/packages/ape/index.html). 

Output File Glossary <a id="output"></a>
--------------------

**sequence\_files/combined\_fasta/combined\_seqs.fna**: The FASTA file containing all sequences, with FASTA headers in the format "&gt;SAMPLEID\_SEQUENCENUMBER".

**nwcs\_arct\_metadata.tsv**: Tab-separated table file containing sample IDs as rows, and metadata as columns.

**reference\_data/97\_otu\_taxonomy.txt**: GreenGenes 13\_8 taxonomic classification file.

**reference\_data/97\_otus.bursttrie\_0.dat, reference\_data/97\_otus.kmer\_0.dat, reference\_data/97\_otus.pos\_0.dat, reference\_data/97\_otus.stats**: GreenGenes 13\_8 taxonomic classification file formatted for SortMeRNA using indexdb\_rna run on reference\_data/97\_otu\_taxonomy.txt.

**reference\_data/97\_otus.fa**: GreenGenes 13\_8 sequences (clustered at 97%).

**reference\_data/core\_set\_aligned.fasta.imputed**: Template for aligning 16S sequences.

**sequence\_files/<SampleID>.assembled.fastq**: FASTQ file containing the sequences successfully assembled by PEAR.

**clustering/rep\_set.fa**: Representative sequences for each OTU. The FASTA labels correspond to the OTU number in the OTU table.

**clustering/final\_otu\_map\_mc1.txt**: File that maps sequences to OTU numbers.

**clustering/pynast\_aligned\_seqs/rep\_set\_aligned.fasta**: FASTA file with aligned OTU representative sequences.

**clustering/pynast\_aligned\_seqs/rep\_set\_aligned\_pfiltered.fasta**: Aligned OTU sequences with empty columns removed.

**clustering/rep\_set.tre**: Phylogenetic tree based on OTU sequences.

**clustering/uclust\_assigned\_taxonomy/rep\_set\_tax\_assignments.txt**: List of OTUs to their taxonomic classification.

**clustering/otu\_table\_high\_conf.tsv**: Tab-delimited OTU table format, low abundance OTUs removed.

**clustering/otu\_table\_high\_conf.biom**: BIOM format (non-human readable) OTU table, low abundance OTUs removed.

**otu\_table/otu\_table\_high\_conf\_rarefied.tsv**: Tab-delimited rarefied OTU table.

**otu\_table/otu\_table\_high\_conf\_rarefied.biom**: BIOM format (non-human readable) rarefied OTU table.

**qiime\_pcoa/pcoa\_plots/pcoa\_2D\_PCoA\_plots.html**: QIIME Bray-Curtis PCoA plot.

Assignment 1 Questions <a id="questions"></a>
----------------------

1) Paired-end assembly: PEAR uses a statistical test to generate a p-value for each paired-end sequence assembly. If the p-value is larger than the cutoff, the sequence will not be assembled. What is the default p-value cutoff?

2) Sequence clustering: We are using the QIIME open-reference pipeline to create clusters of our 16S rRNA sequences. According to [the documentation](http://qiime.org/scripts/pick_open_reference_otus.html), its steps are: (A) assign reads to existing OTUs (closed-reference), (B) take a percentage of the unassigned reads and cluster them de novo, (C) compare all the unassigned reads from step A against the representative sequences from the OTUs generated in step B and assign matches to OTUs, (D) cluster any remaining unassigned reads de novo. Take a look at the results generated in the `clustering` folder. Does it look like all these steps were performed? Is there a particular file that will tell you exactly what was done?

3) Sequence clustering: How many reads were clustered de novo?

3) Taxonomic classification: The ribosomal database project (RDP) is the name for both a taxonomic classifier and a reference dataset. In this lab we have used GreenGenes 13\_8 revision as the reference dataset, and the naive Bayes RDP classifier. What other datasets and classifiers are available? For what reasons might one choose one of these reference datasets?

4) Phylogenetic tree generation: From the clustering log file, can you determine the method used to create a multiple sequence alignment of the OTU sequences? What about if you check the documentation for the QIIME script used to align sequences? What type of multiple sequence alignment algorithm is this? What program was used for the creation of a phylogenetic tree? Which version of the tree building program did we use, and what is the most recent version?

5) OTU table: In this pipeline, QIIME combined a list that maps sequences to OTU names and the taxonomic classifications and create an OTU table. By default, QIIME creates OTU tables in BIOM format (http://biom-format.org/). This is a binary storage format that stores the BIOM file in a sparse representation (no zeros are stored, only positive counts). This results in a smaller file size. We have used the `biom convert` command to create tab-separated OTU tables. When might this less efficient format be preferable?

6) Rarefaction: To eliminate effects due to differing sample sizes, we often subsample the columns of an OTU table without replacement until each column contains the same number of sequences. This process is known more commonly as rarefaction. There has been much debate about the validity of rarefaction. Why might rarefaction be a bad idea? What is an alternate procedure to aid in comparing unevenly sequenced samples?

7) Ordinations: We have created low-dimensional visualizations of our data based on Bray-Curtis and UniFrac distances. PCoA plots attempted to capture as much of the variance present in the distance matrix as is possible in two dimensions. The amount of variance captured can be calculated from the eigenvalues of the covariance matrix calculated from the distance matrix. What percent of the variance was explained by the 3D weighted unifrac PCoA plot? Use this address: http://cbwxx.dyndns.info with xx replaced with your student number to navigate to the results of the analysis (qiime_pcoa_3D/weighted_unifrac_emperor_pcoa_plot/). Can you find any interesting patterns in any of the ordinations? 

8) Exploring the data: Some questions require a bit of mucking through data files to be answered. Being able to navigate through the maze of files you have created is essential. In your OTU table (either the rarefied or the normal OTU table), there should be a relatively large OTU of [genus *Sediminicola*](http://www.bacterio.net/sediminicola.html). Find this OTU in your OTU table, and then track down the sequence of this OTU. Using [NCBI's BLAST tool](http://blast.ncbi.nlm.nih.gov/Blast.cgi), align this sequence against the NR database and exclude uncultured organisms. Does the classification given by RDP seem correct?

[Integrated assignment Day 1 Answer Key](http://bioinformatics-ca.github.io/analysis_of_metagenomic_data_integrated_assignment1_2016/)
------------------------------------------------------------------------------------------------

Integrated Assignment - Day \#2 <a id="day_2"></a>
===============================

In today's session we will try to analyze the metagenomic samples to answer the questions "Who is out there?" and "What are they doing?". We begin by setting up our workspace folder called "integrated\_assignment\_day2" to hold all our analysis files for today.

-   Check if there is already a directory called "integrated\_assignment\_day2"

```
ls -ltrh
rm -rf assignment2
```

-   Make a new folder

```
mkdir -p ~/workspace/assignment2
```

-   Change directory to the new folder

```
cd ~/workspace/assignment2
```

-   To avoid copying all the sample files from the ~/CourseData directory, we will just create a link to that folder from our folder that we just created.

```
ln -s ~/CourseData/integrated_assignment_day2/* .
```

-   Check to see if we have got all the files

```
ls -ltrh 
```

Taxonomic composition of communities in the OSD samples OR "Who is out there?" <a id="who"></a>
------------------------------------------------------------------------------

We will use the MetaPhlan2 program for identifying and quantifying the microbial taxa in our samples.

### Running metaphlan in batch mode

We have the helper script run\_metaphlan2.pl for running MetaPhlan on all our samples one after the other.

-   Run the helper script run\_metaphlan2.pl

```
run_metaphlan2.pl -p 4 -o osd_metaphlan_merged_all.txt *.fasta
```

This step should take ~20 minutes to complete. So, we will try and dissect the run\_metaphlan2.pl PERL script while it is running by opening up another connection to the cloud.

### Dissecting the run\_metaphlan2.pl PERL script

-   Locate the script and open it for viewing; When you open the script using 'less' you can use the up and down arrow to navigate.

```
which run_metaphlan2.pl
```

-   open it for viewing

```
less /usr/local/microbiome_helper-master/run_metaphlan2.pl
```

-   The following are just parts of the script displayed to explain how it works. These are not meant to be copy pasted and run on the terminal.

Getting the options from the user and storing them in variables

```
my ($final_out_file,$parallel,$help);
my $res = GetOptions("output=s" => \$final_out_file,
                     "location=s"=> \$metaphlan_dir,
                     "align_len=i"=>\$align_len,
                     "parallel:i"=>\$parallel,
                     "bowtie:s"=>\$bowtie,
                     "help"=>\$help,
    )or pod2usage(2);`
pod2usage(-verbose=>2) if $help;`
```

Setting up the default paths

```
my $metaphlan_script=$metaphlan_dir.'metaphlan2.py';
my $metaphlan_db=$metaphlan_dir.'db_v20/mpa_v20_m200';
my $metaphlan_pkl=$metaphlan_dir.'db_v20/mpa_v20_m200.pkl';
my $metaphlan_merge=$metaphlan_dir.'utils/merge_metaphlan_tables.py';
```

Build up and run the metaphlan2 command for each sample one after the other

```
foreach my $name (keys %paired_files)
{
   my $pid = $pm->start and next; 
   my $cat;
   if ($gzipped){
       $cat='zcat';
   }else{
       $cat='cat';
   }
   my $out_file=$metaphlan_out_dir.$name;
   my $cmd=join(' ',$cat,@{$paired_files{$name}});
   $cmd.=" | $metaphlan_script  --input_type $format --mpa_pkl $metaphlan_pkl --bt2_ps $bowtie --min_alignment_len $align_len --bowtie2db $metaphlan_db --no_map > $out_file";
   print $cmd,"\n";
   system($cmd);
   $pm->finish;
}
```

### Examining the Metaphlan2 results

-   Check whether the run was successful by examining the output files

```
ls -ltrh
less osd_metaphlan_merged_all.txt
```

We'll now use the PERL script metaphlan\_to\_stamp.pl to convert the Metaphlan output into input for STAMP

-   Convert the metaphlan output to the input profile format required for STAMP

```
metaphlan_to_stamp.pl osd_metaphlan_merged_all.txt > osd_metaphlan_merged_all.spf
```

The next step in the pipeline is to run the progam Humann to identify and quantify the metabolic processes in the metagenomes. This involves comparing the metagenome sequences to the KEGG database using the program diamond. This is a time consuming step and estimated to take around ~80-90 minutes for all our samples. So before moving onto the statistical analysis of the metaphlan results, it would be advisable to start the "Preparing for running Humann" step. You can then switch to your other login instance for continuing with the pipeline.

### Statistical analysis of the taxa using STAMP

STAMP takes two main files as input the profile data which is a table that contains the abundance of features (i.e. taxonomic or functions) and a group metadata file which provides more information about each of the samples in the profile data file.

The metadata file is called "metadata-file-for-osd-subset-210615.txt" and is located in ~/CourseData/integrated\_assignment\_day2/. Download this file locally to your computer using an appropriate method (e.g. WinSCP, etc.) You will also need the profile data file "osd\_metaphlan\_merged\_all.spf" which we generated fro the metaphlan output in the previous step.

-   Load files in STAMP by going to File-&gt; Load Data; You should load both the profile and the metadata files
-   Change the “Profile level” (top left) to “Genus”, ensure that the Group legend (top right) has been set to “depth”, and that “PCA plot” has been set below the large middle window. You should now be looking at a PCA plot where the samples are colored according to their depths.
-   Now change the group field to “prov\_code” and the PCA will be coloured according to that grouping instead.
-   Now lets test what is significantly different between the groups at the Genus rank. Under the “Multiple groups” dialog on the left, check that ANOVA is being used as the statistical test, and select “None” for the multiple test correction. The box at bottom will say what the “Number of active features” is, using these set of statistics.
-   Explore the different visualizations by changing “PCA plot” to each of the other visualizations. Note that you can change which genera is being visualized by selecting different ones on the right hand side. Also, note that you can check the “Show only active features” to reduce the list to those that are significantly different.
-   You can save any plot image using File -&gt; Save plot
-   Switch to the "Two Groups" dialog and select "White's non parametric t-test" from the "Statistical tests" drop-down; check the “Show only active features”; Repeat the same but this time selecting "Benjamini-Hochberg FDR" from the "Multiple test correction" drop-down

Functional composition of the OSD samples OR "What are they doing?" <a id="what"></a>
-------------------------------------------------------------------

We will be using the program Humann for this purpose.

### Preparing for running Humann

-   Perform BLASTX searches agaisnt the KEGG reference database using the program Diamond (We have a helper script for running this step in batch mode for all of our samples)

```
 run_pre_humann.pl -d ~/CourseData/refs/kegg/kegg.reduced -p 8 -o pre_humann *.fasta
```

-   Check whether the run was successful by examining the output files

```
ls -ltrh pre_humann
less pre_humann/OSD106-0m-depth.comb.qc.masked.dedup.subsample.txt
```

### Running Humann

-   copy the humann program folder (humann-0.99) to your working directory

```
cp -r ~/CourseData/refs/humann-0.99/ ~/workspace/
```

-   Change to the humann directory you just copied

```
cd ~/workspace/humann-0.99
```

-   Move the output files for each sample generated by the pre\_humann step into the "input" folder in the humann-0.99 directory

```
mv ~/workspace/assignment2/pre_humann/*.txt input/
```

-   Run the humann program on all samples

```
scons -j 4
```

A bunch of messages will pass on your screen and it should finish in ~20-30 minutes. All of the output is contained in the “humann-0.99/output” directory. There are MANY output files from human. The ones we care about are called:

01b-hit-keg-cat.txt –&gt; KEGG KOs

04b-hit-keg-mpm-cop-nul-nve-nve.txt -&gt; KEGG Modules

04b-hit-keg-mpt-cop-nul-nve-nve.txt -&gt; KEGG Pathways

These files contain relative abundances for each of these different functional classifications. You can look at the format of these using “less”:

-   To statistically test and visualize these resutls using STAMP we need to convert these files into a format readable by STAMP just like we did with the Metaphlan output. We use the "humann\_to\_stamp.pl" PERL script for this step

```
humann_to_stamp.pl output/04b-hit-keg-mpt-cop-nul-nve-nve.txt > pathways.spf

humann_to_stamp.pl output/04b-hit-keg-mpm-cop-nul-nve-nve.txt > modules.spf

humann_to_stamp.pl output/01b-hit-keg-cat.txt > kos.spf
```

-   Since the sample names in the \*.spf files and in the metadata file do not match (the \*.spf files have the text ".subsample" appendedto the sample name), it mght results in errors when you load them in STAMP. So we have to fix the \*.spf by removing the ".subsample" from them.

```
sed -i 's/\.subsample//g' kos.spf

sed -i 's/\.subsample//g' modules.spf

sed -i 's/\.subsample//g' pathways.spf
```

### Statistical analysis of metabolic differences

-   Load the kos.spf file along with the original metadata-file-for-osd-subset-210615.txt file into STAMP.
-   Compare the Arctic samples to the Northwest Atlantic samples using a Two Group test. Use the default Welch’s t-test with BH FDR. Since the number of features (i.e the KO categories) is very high, we will reduce the p-value cut-off

Day 2 assigment questions <a id="questions2"></a>
-------------------------

-   1) How many total sample files do we have?

<!-- -->

-   2) How many sequences does the sample from OSD station 10 contain?

<!-- -->

-   3) How many samples of each type are there in each of the different Province code categories?

<!-- -->

-   4) In the STAMP analysis of the Metaphlan results, do you see any separation in the samples when the PCA is coloured by Depth?

<!-- -->

-   5) In the STAMP analysis of the Metaphlan results, do you see any separation in the samples when the PCA is coloured by the province codes? If so, describe which PC axis differentiates these samples.

<!-- -->

-   6) In the STAMP analysis of the Metaphlan results, in a “multiple group test” using ANOVA with no multiple test correction how many genera are statistically significant?

<!-- -->

-   7) How many are still significant in the “two group test” using White's non-parametric t-test without and with Benjamini-hochberg FDR for multiple test correction?

<!-- -->

-   8) From the kos.spf file what are the top 3 Modules present in the 1m sample from the Bedford basin (station 152)?

<!-- -->

-   9) In the STAMP analysis of the Humann results (with kos.spf file) using a two group test with no multiple test correction applied how many significant differences are seen between the Arctic and Northwest Atlantic samples?

<!-- -->

-   10) What happens when the p-value cut-off is lowered to 0.01 for Q9?

<!-- -->

-   11) In the STAMP analysis of the Humann results with the kos.spf file, what is the most significantly different KEGG pathway? What is the p-value for this KEGG Pathway?

<!-- -->

-   12) Change the p-value to 0.001 and create an “Extended error bar” plot and save the image as a .png using the File-&gt;Save Plot option.

[Integrated assignment Day 2 Answer Key](http://bioinformatics-ca.github.io/analysis_of_metagenomic_data_integrated_assignment2_2016/)
---------------------------------------------------------------------------------------------------------------------
