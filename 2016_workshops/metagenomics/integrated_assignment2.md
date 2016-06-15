---
layout: post2
permalink: /analysis_of_metagenomic_data_integrated_assignment2_2016/
title: Analysis of Metagenomic Data 2016 Student Page
header1: Analysis of Metagenomic Data 2016
header2: Integrated Assignment 2
image: CBW_Metagenome_icon.jpg
---

Background
==========

This integrated lab assignment session will use a subset of samples collected and sequenced through the MicroB3 consortium (http://www.microb3.eu/osd) as part of the global Ocean Sampling Day event in 2014.

To quote from the OSD website:

"The Ocean Sampling Day (OSD) is a simultaneous sampling campaign of the world’s oceans which took place (for the first time) on the summer solstice (June 21st) in the year 2014. These cumulative samples, related in time, space and environmental parameters, provide insights into fundamental rules describing microbial diversity and function and contribute to the blue economy through the identification of novel, ocean-derived biotechnologies. We see OSD data as a reference data set for generations of experiments to follow in the coming decade."

We will be analyzing samples collected from marine stations from various locations in the North West Atlantic shelf and the Polar Atlantic Arctic province. For the integrated assignment on day 1, we will be using sequences from the V4-V6 region of the 16S rRNA gene from 28 samples (19 from the North West Atlantic shelf and 9 from the Polar Atlantic Arctic). On day 2, we will be using shotgun metagenomic sequences from 25 samples (16 from the North West Atlantic shelf and 9 from the Polar Atlantic Arctic). These samples have been downloaded for you already from the European Bioinformatics Institute metagenomics portal (https://www.ebi.ac.uk/metagenomics/projects/ERP009703) and placed in the directories:

~/CourseData/integrated\_assignment\_day1/

~/CourseData/integrated\_assignment\_day2/

Objective: We have chosen a subset of the OSD samples from 2 regions which are expected to differ from each other in terms of temperature and probably also in terms of anthropogenic impact with the samples from the North-west Atlantic shelf being more impacted by human activity than the Arctic samples. The idea is to investigate if these differences are also reflected in a change in the microbial communities in these samples and if so, what are the microbial taxa or metabolic processes that drive these differences.

Integrated Assignment - Day \#1
===============================

Preamble
--------

After the assignment was run, I sensed there was a desire for a "vanilla" QIIME workflow that didn't use all the nasty BASH tricks. I have put a simpler workflow up on GitHub under "vanilla\_qiime\_workflow.sh": <https://github.com/beiko-lab/CBWMeta2015>

This simpler workflow should showcase how easy QIIME is to use. However, I will repeat the warning that it is often invoking software (in particular FastTree, USEARCH, cd-hit, PyNast, and RDP) that is many years out of date. I encourage you to check out <http://qiime.org/install/install.html> and compare the listed versions with the current releases. You should decide if you can accept running your analysis using old (possibly buggy) software. You should also check out <http://qiime.org/scripts/> and read through the documentation and the (many) parameters available to each script before running it.

Franck also pointed us at a great resource for using UPARSE within QIIME without using so many custom commands. The workflow is available at this location: <http://www.brmicrobiome.org/#!standardsand-protocols/cpbw> (click the Illumina 16S profiling link). Some of the commands need to be tweaked slightly for the newest USEARCH v8, but should work for the most part.

Set up
------

Once you have logged into the server, you will create the workspace for today's assignment.

`mkdir -p ~/workspace/assignment1`
`cd ~/workspace/assignment1`

All of the required commands are in a BASH script which we can copy to our workspaces. Each of the commands contained within this script are displayed and discussed below.

`wget `[`https://raw.githubusercontent.com/beiko-lab/CBWMeta2015/master/Integrated_Lab_1.sh`](https://raw.githubusercontent.com/beiko-lab/CBWMeta2015/master/Integrated_Lab_1.sh)
`chmod u+x Integrated_Lab_1.sh`

The first command fetches the script from the web. The second command allows us to execute it.

This script takes about half an hour to run to completion. When running the script, piping to the \`tee\` command is helpful to store a copy of the output for later viewing.

`./Integrated_Lab_1.sh 2>&1 | tee -a log.txt`

Let's go through the script command by command. First is a bunch of setup commands.

`mkdir sequence_files`
`mkdir reference_data`
`mkdir scripts`
`ln -s ~/CourseData/integrated_assignment_day1/*.fastq.gz sequence_files/`
`ln -s ~/CourseData/integrated_assignment_day1/nwcs_arct_metadata.tsv .`
`ln -s ~/CourseData/integrated_assignment_day1/97* reference_data/`
`ln -s ~/CourseData/integrated_assignment_day1/core_set_aligned.fasta.imputed reference_data/`

These commands get the necessary data files into our workspace. Source FASTQ sequence files are linked to the sequence\_files directory. GreenGenes reference set and alignment template are linked to the reference\_data directory.

`wget -O scripts/mesas-pcoa `[`https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa`](https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-pcoa)
`wget -O scripts/mesas-uc2clust `[`https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-uc2clust`](https://raw.githubusercontent.com/neufeld/MESaS/master/scripts/mesas-uc2clust)
`chmod u+x scripts/*`

Here we download a few scripts we will need from GitHub, and then the \`chmod\` command allows us to execute them.

Workflow
--------

### Paired-end Assembly

The first step in our pipeline is to assemble the Illumina paired-end reads with [PEAR](https://github.com/xflouris/PEAR).

` for i in "1" "5" "9" "13" "17" "21" "25"`
`do`
`    echo $i`
`    find ~/workspace/assignment1/sequence_files -name "*.fastq.gz" -printf '%f\n' | sed 's/_.*//' | sort | uniq | sed -n $i,$((i+3))p | \`
`while read line; do ( pear -f sequence_files/${line}_1.fastq.gz -r sequence_files/${line}_2.fastq.gz -o ${line} & ); done > /dev/null`
`    sleep 60`
`done`

Don't panic! This complex command is meant to showcase the power of the tools available in most BASH shell environments. This command is finding the forward and reverse sequence read files and launching 4 simultaneous PEAR processes every minute. By launching multiple PEAR processes at once, we save processing time. The complicated run-up with \`find\`, \`sed\`, \`sort\`, and \`uniq\` finds the files we want to assemble automatically so that we don't have to type all 54 file names in manually. A simpler approach would be to copy and paste the PEAR command, replacing the forward and reverse file arguments each time.

`while [ $( ps -ef | grep pear | wc -l ) -gt 1 ]`
`do`
` sleep 10`
` echo "Waiting on PEAR to finish..."`
`done`

This command detects if there are any PEAR processes running, and if there are, it will wait for them to finish.

`for filename in $( ls *.assembled.fastq )`
`do`
`    awk 'BEGIN{ORS=""; i=0;}{split(FILENAME, x, "."); prefix=x[1]; sub("@","",prefix); print ">" prefix "_" i "\n"; i+=1; `
`         getline; print; print "\n"; getline; getline;}' ${filename} >> seq.fasta`
`done`

Here we leverage \`awk\` to convert the assembled FASTQ files into FASTA files where the sequence headers are in the format "&gt;SAMPLEID\_SEQUENCENUMBER". This is the format that QIIME expects for its input sequence files.

`mkdir sequence_files`
`mv *.fastq sequence_files`

Staying organized is important! It is easy to create hundreds of files and forget which is which.

### Clustering with UPARSE

Next, we start the [UPARSE pipeline](http://drive5.com/usearch/manual/uparse_pipeline.html), which is available in the USEARCH software. Many clustering workflows are programmed into QIIME and Mothur, but here we are interacting directly with USEARCH to better illustrate the steps involved. At each stage, USEARCH will output useful statistics to the screen. This is why it's always a good idea to save your terminal output!

`usearch -derep_fulllength seq.fasta -fastaout derep.fa -sizeout`

This collapses the sequence file down to only contain unique sequences. This will improve performance in the later stages.

`usearch -sortbysize derep.fa -fastaout sorted.fa -minsize 2`

The sequences are then sorted by abundance (most abundant first), and all sequences that are only observed once (singletons) are temporarily set aside.

`usearch -cluster_otus sorted.fa -otus otus.fa -otu_radius_pct 3 -sizeout -uparseout results.txt`

The abundance-sorted non-singletons are clustered at 97% sequence identity. This results in a FASTA file with the OTU representative sequences, but we don't yet know which sequences belong to which OTUs.

`awk 'BEGIN{count=0;}{if ($0~/>/){print ">" count; count+=1;} else {print}}' otus.fa > rep_set.fasta`

The FASTA labels are modified to be a little more QIIME-friendly. Labels will be &gt;0, &gt;1, &gt;2, ... etc.

`usearch -usearch_global seq.fasta -db rep_set.fasta -strand both -id 0.97 -uc map.uc -threads 4`

Using the USEARCH algorithm, we now map all of our original sequences onto the OTU representative sequences produced in the clustering step.

`scripts/mesas-uc2clust -t 4 map.uc seq_otus.txt`

Finally, we use a custom Python script to convert the USEARCH hit file to the format that QIIME expects for OTU lists (ie, "OTUID\\tSEQ1\\tSEQ2\\t...").

### Classification, Alignment, Trees, and OTU Tables with QIIME

Next, we do some housekeeping to clean up and ensure that we are using the correct software.

`mkdir cluster`
`mv results.txt otus.fa map.uc rep_set.fasta seq_otus.txt sorted.fa derep.fa cluster`
`export PYTHONPATH=~/local/lib/python2.7/site-packages`
`export RDP_JAR_PATH=/usr/local/rdp_classifier_2.2/rdp_classifier-2.2.jar`

The third command ensures we are using the proper BIOM software version. The fourth command points QIIME at the RDP classifier location.

`assign_taxonomy.py -m rdp -i cluster/rep_set.fasta -o taxonomic_classification -t reference_data/97_otu_taxonomy.txt -r reference_data/97_otus.fasta -c 0.6`

This QIIME script will classify our OTU representative sequences (rep\_set.fasta) with RDP against the GreenGenes 13\_8 reference set. Posterior probability of at least 60% is required for a classification to be saved.

`align_seqs.py -m pynast -i cluster/rep_set.fasta -o alignment -t reference_data/core_set_aligned.fasta.imputed`
`filter_alignment.py -i alignment/rep_set_aligned.fasta -o alignment -s`
`mkdir phylogeny`
`make_phylogeny.py -i alignment/rep_set_aligned_pfiltered.fasta -t fasttree -o phylogeny/rep_set.tre -l phylogeny/log.txt`

These QIIME commands take the OTU representative sequence set in, create a multiple sequence alignment, filter empty columns out of it, and finally create a phylogenetic tree. The core\_set\_aligned.fasta.imputed file is used to align our sequences against. It can be downloaded from [this location at the GreenGenes FTP site](ftp://greengenes.microbio.me/greengenes_release/unversioned/).

`mkdir otu_table`
`make_otu_table.py -i cluster/seq_otus.txt -o otu_table/otu_table.biom -t taxonomic_classification/rep_set_tax_assignments.txt`

We provide our cluster list file (seq\_otus.txt) and taxonomic classification file, and QIIME will produce a BIOM formatted OTU table.

`biom convert --to-tsv -i otu_table/otu_table.biom --table-type='OTU table' -o otu_table/otu_table.tab`
` --header-key=taxonomy --output-metadata-id=Consensus\ Lineage`

Using the [biom-format software package](https://github.com/biocore/biom-format), we convert the BIOM table to a more traditional format.

`single_rarefaction.py -i otu_table/otu_table.biom -o otu_table/otu_table_rarefied.biom -d `
``  `awk 'BEGIN{FS="\t"} NR == 1 { } NR == 2 { max = NF-1; } NR > 2 { for (i = 2; i <= max; i++) { c[i] += $i; } } ``
``  END { smallest = c[2]; for (i = 3; i <= max; i++) { if (c[i] < smallest) { smallest = c[i]; }} print smallest; }' otu_table/otu_table.tab` ``

This complex command consists of two parts. The embedded \`awk\` command will traverse our newly converted tab-separated OTU table to find the minimum number of sequences in a sample. This is passed as an argument to QIIME's rarefaction script, which produces a BIOM format rarefied OTU table.

`biom convert --to-tsv -i otu_table/otu_table_rarefied.biom --table-type='OTU table' -o otu_table/otu_table_rarefied.tab`
` --header-key=taxonomy --output-metadata-id=Consensus\ Lineage`

Like before, we can convert the rarefied OTU table from BIOM format to tab-separated format.

### Downstream Analyses

We have our OTU tables! All our hard work has paid off, and we can now explore the data with visualizations and other analyses.

`#summarize_taxa_through_plots.py -f -s -i otu_table/otu_table.biom -o taxaplot -m nwcs_arct_metadata.tsv`

This command is commented out because it tends to take a long time to run. It will create taxonomy bar plots for each sample at each level from kingdom to genus. If you wish to run this after the script has completed, you will first need to run \`export PYTHONPATH=~/local/lib/python2.7/site-packages\` to point your session at the proper BIOM format software.

Now we create the every popular PCoA plots.

`beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m weighted_unifrac -t phylogeny/rep_set.tre`
`beta_diversity.py -i otu_table/otu_table_rarefied.biom -o qiime_pcoa/distance/ -m bray_curtis`

The first step is to create the distance matrices. Note that for UniFrac, a tree is required since it is a phylogenetic beta diversity measure.

`mkdir qiime_pcoa/pcoa`
`principal_coordinates.py -i qiime_pcoa/distance/ -o qiime_pcoa/pcoa`

Next, we run the PCoA algorithm on the data matrices. This will rotate and scale the data in order to capture as much variance as possible in two (or three) dimensions.

`mkdir qiime_pcoa/pcoa_plots`
`make_2d_plots.py -i qiime_pcoa/pcoa -m nwcs_arct_metadata.tsv -o qiime_pcoa/pcoa_plots/`

QIIME will also plot this for us, and create a unique colouring for each metadata category in the nwcs\_arct\_metadata.tsv mapping file.

`scripts/mesas-pcoa -i otu_table/otu_table_rarefied.tab -m nwcs_arct_metadata.tsv -d bray -o nmds`

Finally, we try the same PCoA with a different method: an Rscript. mesas-pcoa is an R script that uses the [R vegan package](http://cran.r-project.org/web/packages/vegan/index.html) and [ape package](http://cran.r-project.org/web/packages/ape/index.html). These packages contain many different functions useful for ecology.

Output File Glossary
--------------------

**seq.fasta**: The FASTA file containing all sequences, with FASTA headers in the format "&gt;SAMPLEID\_SEQUENCENUMBER".

**nwcs\_arct\_metadata.tsv**: Tab-separated table file containing sample IDs as rows, and metadata as columns.

**reference\_data/97\_otu\_taxonomy.txt**: GreenGenes 13\_8 taxonomic classification file.

**reference\_data/97\_otus.fa**: GreenGenes 13\_8 sequences (clustered at 97%).

**reference\_data/core\_set\_aligned.fasta.imputed**: Template for aligning 16S sequences.

**sequence\_files/<SampleID>.assembled.fastq**: FASTQ file containing the sequences successfully assembled by PEAR.

**cluster/derep.fa**: Contains only unique sequences.

**cluster/rep\_set.fa**: Representative sequences for each OTU. The FASTA labels correspond to the OTU number in the OTU table.

**cluster/seq\_otus.txt**: File that maps sequences to OTU numbers.

**alignment/rep\_set\_aligned.fasta**: FASTA file with aligned OTU representative sequences.

**alignment/rep\_set\_aligned\_pfiltered.fasta**: Aligned OTU sequences with empty columns removed.

**phylogeny/rep\_set.tre**: Phylogenetic tree based on OTU sequences.

**taxonomic\_classification/rep\_set\_tax\_assignments.txt**: List of OTUs to their taxonomic classification.

**otu\_table/otu\_table.tab**: Tab-delimited OTU table format.

**otu\_table/otu\_table.biom**: BIOM format (non-human readable) OTU table.

**otu\_table/otu\_table\_rarefied.tab**: Tab-delimited rarefied OTU table.

**otu\_table/otu\_table\_rarefied.biom**: BIOM format (non-human readable) rarefied OTU table.

**qiime\_pcoa/pcoa\_plots/pcoa\_2D\_PCoA\_plots.html**: QIIME Bray-Curtis PCoA plot.

**mesas\_pcoa/pcoa-bray.pdf**: PCoA generated by R packages.

**mesas\_pcoa/eigenvalues.txt**: Eigenvalues of PCoA generated by R packages.

Assignment 1 Questions
----------------------

1) Paired-end assembly: PEAR uses a statistical test to generate a p-value for each paired-end sequence assembly. If the p-value is larger than the cutoff, the sequence will not be assembled. What is the default p-value cutoff?

2) Sequence clustering: We are using the UPARSE pipeline to create de novo clusters of our 16S rRNA sequences. Its steps are: collapse unique sequences, sort by abundance, cluster all unique sequences seen more than once while checking for chimeras, then map all sequences back onto the OTUs. How many chimeras did the UPARSE pipeline detect? Approximately how many sequences were discarded by the UPARSE pipeline? How many OTUs were generated?

3) Taxonomic classification: The ribosomal database project (RDP) is the name for both a taxonomic classifier and a reference dataset. In this lab we have used GreenGenes 13\_8 revision as the reference dataset, and the naive Bayes RDP classifier. What other datasets and classifiers are available? For what reasons might one choose one of these reference datasets?

4) Phylogenetic tree generation: From the script file, can you determine the method used to create a multiple sequence alignment of the OTU sequences? What type of multiple sequence alignment algorithm is this? What program was used for the creation of a phylogenetic tree? Which version of the tree building program did we use, and what is the most recent version?

5) OTU table: In this pipeline, QIIME combined a list that maps sequences to OTU names and the taxonomic classifications and create an OTU table. By default, QIIME creates OTU tables in BIOM format (http://biom-format.org/). This is a binary storage format that stores the BIOM file in a sparse representation (no zeros are stored, only positive counts). This results in a smaller file size. We have used the \`biom convert\` command to create tab-separated OTU tables. When might this less efficient format be preferable?

6) Rarefaction: To eliminate effects due to differing sample sizes, we often subsample the columns of an OTU table without replacement until each column contains the same number of sequences. This process is known more commonly as rarefaction. There has been much debate about the validity of rarefaction. Why might rarefaction be a bad idea? What is an alternate procedure to aid in comparing unevenly sequenced samples?

7) Ordinations: We have created low-dimensional visualizations of our data based on Bray-Curtis and UniFrac distances. PCoA plots attempted to capture as much of the variance present in the distance matrix as is possible in two dimensions. The amount of variance captured can be calculated from the eigenvalues of the covariance matrix calculated from the distance matrix. What percent of the variance was explained by the Bray-Curtis PCoA plot? Can you find any interesting patterns in any of the ordinations? **Note:** This question had previously reference an NMDS plot that was not included in the workflow. Ignore that, sorry!

8) Exploring the data: Some questions require a bit of mucking through data files to be answered. Being able to navigate through the maze of files you have created is essential. In your OTU table (either the rarefied or the normal OTU table), there should be a relatively large OTU of [genus *Sediminicola*](http://www.bacterio.net/sediminicola.html). Find this OTU in your OTU table, and then track down the sequence of this OTU. Using [NCBI's BLAST tool](http://blast.ncbi.nlm.nih.gov/Blast.cgi), align this sequence against the NR database and exclude uncultured organisms. Does the classification given by RDP seem correct?

[Integrated assignment Day 1 Answer Key](Metagenomics_IntegratedAssignment_AnswerKey "wikilink")
------------------------------------------------------------------------------------------------

Integrated Assignment - Day \#2
===============================

In today's session we will try to analyze the metagenomic samples to answer the questions "Who is out there?" and "What are they doing?". We begin by setting up our workspace folder called "integrated\_assignment\_day2" to hold all our analysis files for today.

-   Check if there is already a directory called "integrated\_assignment\_day2"

`ls -ltrh`
`rm -rf assignment2`

-   Make a new folder

`mkdir -p ~/workspace/assignment2`

-   Change directory to the new folder

`cd ~/workspace/assignment2`

-   To avoid copying all the sample files from the ~/CourseData directory, we will just create a link to that folder from our folder that we just created.

`ln -s ~/CourseData/integrated_assignment_day2/* .`

-   Check to see if we have got all the files

`ls -ltrh `

Taxonomic composition of communities in the OSD samples OR "Who is out there?"
------------------------------------------------------------------------------

We will use the MetaPhlan2 program for identifying and quantifying the microbial taxa in our samples.

### Running metaphlan in batch mode

We have the helper script run\_metaphlan2.pl for running MetaPhlan on all our samples one after the other.

-   Run the helper script run\_metaphlan2.pl

`run_metaphlan2.pl -p 4 -o osd_metaphlan_merged_all.txt *.fasta`

This step should take ~20 minutes to complete. So, we will try and dissect the run\_metaphlan2.pl PERL script while it is running by opening up another connection to the cloud.

### Dissecting the run\_metaphlan2.pl PERL script

-   Locate the script and open it for viewing; When you open the script using 'less' you can use the up and down arrow to navigate.

`which run_metaphlan2.pl`

-   open it for viewing

`less /usr/local/microbiome_helper-master/run_metaphlan2.pl`

-   The following are just parts of the script displayed to explain how it works. These are not meant to be copy pasted and run on the terminal.

Getting the options from the user and storing them in variables

`my ($final_out_file,$parallel,$help);`
`my $res = GetOptions("output=s" => \$final_out_file,`
`                     "location=s"=> \$metaphlan_dir,`
`                     "align_len=i"=>\$align_len,`
`                     "parallel:i"=>\$parallel,`
`                     "bowtie:s"=>\$bowtie,`
`                     "help"=>\$help,`
`    )or pod2usage(2);`
`pod2usage(-verbose=>2) if $help;`

Setting up the default paths

`my $metaphlan_script=$metaphlan_dir.'metaphlan2.py';`
`my $metaphlan_db=$metaphlan_dir.'db_v20/mpa_v20_m200';`
`my $metaphlan_pkl=$metaphlan_dir.'db_v20/mpa_v20_m200.pkl';`
`my $metaphlan_merge=$metaphlan_dir.'utils/merge_metaphlan_tables.py';`

Build up and run the metaphlan2 command for each sample one after the other

`foreach my $name (keys %paired_files)`
`{`
`   my $pid = $pm->start and next; `
`   my $cat;`
`   if ($gzipped){`
`       $cat='zcat';`
`   }else{`
`       $cat='cat';`
`   }`
`   my $out_file=$metaphlan_out_dir.$name;`
`   my $cmd=join(' ',$cat,@{$paired_files{$name}});`
`   $cmd.=" | $metaphlan_script  --input_type $format --mpa_pkl $metaphlan_pkl --bt2_ps $bowtie --min_alignment_len $align_len --bowtie2db $metaphlan_db --no_map > $out_file";`
`   print $cmd,"\n";`
`   system($cmd);`
`   $pm->finish;`
`}`

### Examining the Metaphlan2 results

-   Check whether the run was successful by examining the output files

`ls -ltrh`
`less osd_metaphlan_merged_all.txt`

We'll now use the PERL script metaphlan\_to\_stamp.pl to convert the Metaphlan output into input for STAMP

-   Convert the metaphlan output to the input profile format required for STAMP

`metaphlan_to_stamp.pl osd_metaphlan_merged_all.txt > osd_metaphlan_merged_all.spf`

The next step in the pipeline is to run the progam Humann to identify and quantify the metabolic processes in the metagenoems. This involves comparing the metagenome sequences to the KEGG database using the prgram diamond. This is a time comsuming step and estimated to take around ~80-90 minutes for all our samples. So before moving onto the statistical analysis of the metaphlan results, it would be advisable to start the "Preparing for running Humann" step. You can then switch to your other login instance for continuing with the pipeline.

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

Functional composition of the OSD samples OR "What are they doing?"
-------------------------------------------------------------------

We will be using the program Humann for this purpose.

### Preparing for running Humann

-   Perform BLASTX searches agaisnt the KEGG reference database using the program Diamond (We have a helper script for running this step in batch mode for all of our samples)

` run_pre_humann.pl -d ~/CourseData/refs/kegg/kegg.reduced -p 8 -o pre_humann *.fasta`

-   Check whether the run was successful by examining the output files

`ls -ltrh pre_humann`
`less pre_humann/OSD106-0m-depth.comb.qc.masked.dedup.subsample.txt`

### Running Humann

-   copy the humann program folder (humann-0.99) to your working directory

`cp -r ~/CourseData/refs/humann-0.99/ ~/workspace/`

-   Change to the humann directory you just copied

`cd ~/workspace/humann-0.99`

-   Move the output files for each sample generated by the pre\_humann step into the "input" folder in the humann-0.99 directory

`mv ~/workspace/assignment2/pre_humann/*.txt input/`

-   Run the humann program on all samples

`scons -j 4`

A bunch of messages will pass on your screen and it should finish in ~20-30 minutes. All of the output is contained in the “humann-0.99/output” directory. There are MANY output files from human. The ones we care about are called:

01b-hit-keg-cat.txt –&gt; KEGG KOs

04b-hit-keg-mpm-cop-nul-nve-nve.txt -&gt; KEGG Modules

04b-hit-keg-mpt-cop-nul-nve-nve.txt -&gt; KEGG Pathways

These files contain relative abundances for each of these different functional classifications. You can look at the format of these using “less”:

-   To statistically test and visualize these resutls using STAMP we need to convert these files into a format readable by STAMP just like we did with the Metaphlan output. We use the "humann\_to\_stamp.pl" PERL script for this step

`humann_to_stamp.pl output/04b-hit-keg-mpt-cop-nul-nve-nve.txt > pathways.spf`

`humann_to_stamp.pl output/04b-hit-keg-mpm-cop-nul-nve-nve.txt > modules.spf`

`humann_to_stamp.pl output/01b-hit-keg-cat.txt > kos.spf`

-   Since the sample names in the \*.spf files and in the metadata file do not match (the \*.spf files have the text ".subsample" appendedto the sample name), it mght results in errors when you load them in STAMP. So we have to fix the \*.spf by removing the ".subsample" from them.

`sed -i 's/\.subsample//g' kos.spf`

`sed -i 's/\.subsample//g' modules.spf`

`sed -i 's/\.subsample//g' pathways.spf`

### Statistical analysis of metabolic differences

-   Load the kos.spf file along with the original metadata-file-for-osd-subset-210615.txt file into STAMP.
-   Compare the Arctic samples to the Northwest Atlantic samples using a Two Group test. Use the default Welch’s t-test with BH FDR. Since the number of features (i.e the KO categories) is very high, we will reduce the p-value cut-off

Day 2 assigment questions
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

![](pdf.gif "fig:pdf.gif") [Integrated assignment Day 2 Answer Key](Media:Cbw-workshop-assignment-key.pdf "wikilink")
---------------------------------------------------------------------------------------------------------------------
