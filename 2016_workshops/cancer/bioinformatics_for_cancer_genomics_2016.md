---
layout: post2
permalink: /bioinformatics_for_cancer_genomics_2016/
title: Bioinformatics for Cancer Genomics 2016 Student Page
header1: Bioinformatics for Cancer Genomics 2016
header2: Workshop pages for students
image: CBW_cancerDNA_icon-16.jpg
---

#### Contents
[Course Schedule](#course_schedule)

[Workshop Q/A Forum](#q_a_forum)

[Laptop Setup Instructions](#laptop_setup)

[Pre-Workshop Tutorials](#pre_tutorials)

[Pre-Workshop Readings](#pre_readings)

[Logging into the Amazon Cloud](#amazon_cloud)

...[Logging in with ssh (Mac/Linux)](#ssh_login)

...[Logging in with Putty (Windows)](#putty_login)

...[File System Layout](#file_system_layout)

**[Day 1](#day_1)**


  ...[Welcome](#welcome)
  
  ...[Module 1: Introduction to Cancer Genomics](#module_1)
  
  ...[Module 2.1: Databases and Visualization Tools](#module_2.1)
  
  ...[Module 2.1: Logging into the Cloud](#module_2.2)
  
  ...[R Review Session](#r_review)


**[Day 2](#day_2)**


  ...[Module 3: Mapping and Genome Rearrangement](#module_3)
  
  ...[Module 4: Gene Fusion Discovery](#module_4)
  

**[Day 3](#day_3)**

  
  ...[Module 5: Copy Number Alterations](#module_5)
  
  ...[Module 6: Somatic Mutations](#module_6)
  

**[Day 4](#day_4)**


  ...[Module 7: Gene Expression Profiling](#module_7)
  
  ...[Module 8: Variants to Networks](#module_8)
  
  ......[Part 1: How to annotate variants and prioritize potentially relevant ones ](#part_1)
  
  ......[Part 2: How to annotate variants and prioritize potentially relevant ones ](#part_2)
  
  
  **[Day 5](#day_5)**
  
  
  ......[Part 3: How to annotate variants and prioritize potentially relevant ones ](#part_3)
  
  ...[Module 9: Integration of Clinical Data ](#module_9)
  
  

***
***

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/BiCG_2016_Schedule_v1.pdf">Schedule for May 30 to June 3, 2016</a>


###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/CancerGenomics2016">here</a>!


###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/laptop_setup_instructions.pdf">here</a>.
  
##### Difference Between **R** and **RStudio**
<br>
**RStudio** doesn't know where libraries are installed, when they are not installed through the **RStudio** package manager. To tell **RStudio** the location, you can define the path in a startup file. Create a file called .Renviron . Inside there:

```r
R_LIBS=<R Library Path of other installed packages>
```

That was the problem when students installed things in **RStudio** at the command line using the **R** command <code>install.package()</code>.

... or you could use the package manger to install libraries.

##### Syntax highlighting
<br>
... of scripts in the R editor does not seem to work under Windows. If you want highlighted syntax, use RStudio instead.


###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [CBW R tutorial](http://bioinformatics.ca/workshop_wiki/index.php/R_tutorial) or [R Tutorial](http://www.cyclismo.org/tutorial/R/) 
* The [R command cheat sheet](../../resources/R_Short-refcard.pdf)
* [PlottingReference.pdf](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/resources/Plotting.Reference.ipynb)


2) **Cytoscape 3.x Preparation tutorials**: Complete the [introductory tutorial to Cytoscape 3.x](http://opentutorials.cgl.ucsf.edu/index.php/Portal:Cytoscape3): 
* Introduction to Cytoscape3 - User Interface
* Introduction to Cytoscape3 - Welcome Screen
* Introduction to Cytoscape 3.1 - Networks, Data, Styles, Layouts and App Manager


3) **UNIX Preparation tutorials**: Please complete tutorials #1-3 on [UNIX Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)
* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 



###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop, read these.
  
[Database resources of the National Center for Biotechnology Information](http://www.ncbi.nlm.nih.gov/pubmed/26615191)

[COSMIC: mining complete cancer genomes in the Catalogue of Somatic Mutations in Cancer](http://www.ncbi.nlm.nih.gov/pubmed/20952405/)

[Integrative genomic profiling of human prostate cancer](http://www.ncbi.nlm.nih.gov/pubmed/20579941)

[Predicting the functional impact of protein mutations: application to cancer genomics](http://www.ncbi.nlm.nih.gov/pubmed/21727090)

[Cancer genome sequencing study design](http://www.ncbi.nlm.nih.gov/pubmed/23594910)

[Using cloud computing infrastructure with CloudBioLinux, CloudMan, and Galaxy](http://www.ncbi.nlm.nih.gov/pubmed/22700313)

[The UCSC Genome Browser database: extensions and updates 2013](http://www.ncbi.nlm.nih.gov/pubmed/23155063)

[Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration](http://www.ncbi.nlm.nih.gov/pubmed/22517427)

[Feature-based classifiers for somatic mutation detection in tumour–normal paired sequencing data](http://www.ncbi.nlm.nih.gov/pubmed/22084253)

[Expression Data Analysis with Reactome](http://www.ncbi.nlm.nih.gov/pubmed/25754994)

  
  
### Logging into the Amazon Cloud <a id="amazon_cloud"></a>

* These instructions will **ONLY** be relevant in class, as the Cloud will not be accessible from home in advance of the class.

* We have set up 30 instances on the Amazon cloud - one for each student. In order to log in to your instance, you will need a security certificate. If you plan on using Linux or Mac OS X, please download this **link certificate here**. Otherwise if you plan on using Windows (with Putty and Winscp), please download this **link certificate here**. 

* On the cloud, we're going to use the default username: **ubuntu**

#### Logging in with ssh (Mac/Linux) <a id="ssh_login"></a>
<p>
<h5> Logging in </h5>
</p>

* Make sure the permissions on your certificate are secure. Use chmod on your downloaded certificate:

```bash
 chmod 600 CBWCG.pem
```

* To log in to the node, use the -i command line argument to specify your certificate:

```bash
 ssh -i CBWCG.pem ubuntu@cbw#.dyndns.info
```

(where # is your assigned student number. Your student number is the number on the participant list. If your number is less than 10, please add 0 in front of it.)

##### Copying files to your computer
<br>
* To copy files from an instance, use scp in a similar fashion:

```bash
 scp -i CBWCG.pem ubuntu@cbw#.dyndns.info:CourseData/genome/g1k/human_g1k_v37.fasta.fai .
```

* Everything created in your workspace on the cloud is also available by a web server on your cloud instance.  Simply go to the following in your browser:

 http://cbw#.dyndns.info/ http://cbw#.dyndns.info/

#### Logging in with Putty (Windows) <a id="putty_login"></a>
<p>
<h5> Logging in </h5>
</p>
To configure Putty, start Putty and do the following:

* Fill in the "Host name" field with cbw#.dyndns.info (where # is your assigned student number. Your student number is the number on the participant list. If your number less is than 10, please add 0 in front of it.)

* In the left hand categories,under the Connection category choose Data.  In the auto-login username field write ***ubuntu***.

* In the left hand categories, in the Connection category next to SSH click on the **+**. Click on Auth. In the private-key file for authentication field, hit browse and find the CBWCG.ppk certificate that you downloaded above.

* In the left hand categories, click on Session.  In the Saved Sessions field write **Amazon node** and click save.

**Now that Putty is configured**, all you have to do is start putty and double-click on "Amazon node" to login.

##### Copying files to your computer
<br>
To configure WinScp, start WinScp and do the following:

* On the right-hand buttons click "New".

* Fill in the "Host name" field with cbw#.dyndns.info (where # is your assigned student number. Your student number is the number on the participant list. If your number is less than 10, please add 0 in front of it.)

* Fill in the "User name" field with **ubuntu**

* Leave the password field empty

* In the "private key file" field press the "..." button to browse for the CBWCG.ppk certificate.

* Click the "Save..." button and in the "Save session as" field write "Amazon node" .


**Now that WinScp is configured**, all you have to do is start WinScp and double-click on **Amazon node** to start copying files.

#### File System Layout <a id="file_system_layout"></a>
<br>
When you log in, you'll notice that you have two directories: **CourseData** and **workspace**.

* The **CourseData** directory will contain the files that you'll need to complete your lab assignments.

* The **workspace** directory is where we will keep our temporary files. By default, we have around 40 GB available for our output files in the workspace directory. If you run out of space, you may need to delete some files from this directory.




***

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="green">Ann Meyer</font>* 
<br>

***

###  Module 1: Introduction to Cancer Genomics <a id="module_1"></a>

  *<font color="green">TBA</font>*
  
  Lecture:


***

###  Module 2.1: Databases and Visualization Tools <a id="module_2.1"></a>

  *<font color="green">Francis Ouellette</font>*
  
  Lecture:
  
  Lab practical:
  
#### Links:
  
 * [ICGC](http://icgc.org/) 
 * [DCC portal on ICGC](http://dcc.icgc.org/) 
 * [Docs for ICGC](http://docs.icgc.org/) 
 * [Integrated Genomics Viewer](http://www.broadinstitute.org/igv/) 
 * [UCSC Genome Browser](http://genome.ucsc.edu/) 
 * [UCSC Genome Browser](https://genome-cancer.ucsc.edu/) 
 * [Cancer Genome Workbench](https://cgwb.nci.nih.gov/) 
 * [cBioPortal for Cancer Genomics](http://www.cbioportal.org/public-portal/web_api.jsp/) 
 * [Savant Genome Browser](http://genomesavant.com/p/home/index/)
  

***  

###  Module 2.2: Logging into the Cloud <a id="module_2.2"></a>

  *<font color="green">Francis Ouellette</font>*
  
  Lecture:


***

###  *Optional* **R** Review Session <a id="r_review"></a>

*<font color="green">TBA</font>*

Lecture:

Lab practical:

#### Links:

 * [R Statistical Package](http://www.r-project.org/) 
 * [R Studio](http://rstudio.org/) 


***

##  Day 2 <a id="day_2"></a>

###  Module 3: Mapping and Genome Rearrangement <a id="module_3"></a>

  *<font color="green">Jared Simpson</font>*
  
  Lecture:
  
  Lab practical:
  
  **Bonus**
You can view your results (BAM and BAM.BAI file) in the IGV browser by using the URL for that file from your Cloud instance. We have a web server running on the Amazon cloud for each instance. In a browser, like Firefox, type in your server name (cbw#.dyndns.info) and all files under your workspace will be shown there. Find your Bam and Bam.Bai file, right click it and 'copy the link location'. Start IGV and choose 'load from URL' from File menu, and then paste the location you just copied and you will see the Bam file you just generated in IGV! Narrow down the view to chromosome 15 or 17 where the break points were identified.

#### Links:

 * [What does my SAM flag mean?](https://broadinstitute.github.io/picard/explain-flags.html)
 * [Tools for Mapping High-throughput Sequencing Data Paper (2012)](http://www.ncbi.nlm.nih.gov/pubmed/23060614) 
 * [SAM/BAM file specifications](http://samtools.sourceforge.net/SAM1.pdf) 
 * [samtools](http://samtools.sourceforge.net/) 
 * [Picard](http://picard.sourceforge.net) 
 * [bwa](http://bio-bwa.sourceforge.net/) 
 * [GASV](http://compbio.cs.brown.edu/software.html) 
 * [BreakDancer](http://breakdancer.sourceforge.net/) 
  
##### Extras:

* exampleheader.txt **Put example header here**
* [explain sam flag app](https://broadinstitute.github.io/picard/explain-flags.html) 
* [sam flags explained](http://davetang.org/muse/2014/03/06/understanding-bam-flags/) 

***

###  Module 4: Gene Fusion Discovery <a id="module_4"></a>

  *<font color="green">Andrew McPherson</font>*
  
  Lecture:
  
  Lab practical:
  
#### Papers and Background Material:


 * [Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.ncbi.nlm.nih.gov/pubmed/22383036)
 
 * [ENCODE RNA-seq Standards](https://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf)


#### Links:

 * [BioStar](http://www.biostars.org/)
 * [SeqAnswers](http://seqanswers.com/)
 * [Integrative Genomics Viewer (IGV)](http://www.broadinstitute.org/igv/)
 * [FASTQ format](http://en.wikipedia.org/wiki/FASTQ_format)
 * [SAM/BAM format](http://samtools.sourceforge.net/SAM1.pdf)
 * [Illumina iGenomes](http://tophat.cbcb.umd.edu/igenomes.html)
 * [SamTools](http://samtools.sourceforge.net/)
 * [Picard](http://picard.sourceforge.net/)
 * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
 * [SAMStat](http://samstat.sourceforge.net/)
 * [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
 * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
 * [TopHat/TopHat2](http://tophat.cbcb.umd.edu/)
 * [Cufflinks/Cuffdiff](http://cufflinks.cbcb.umd.edu/)
 * [CummeRbund](http://compbio.mit.edu/cummeRbund/)


***

##  Day 3 <a id="day_3"></a>

###  Module 5: Copy Number Alterations <a id="module_5"></a>

  *<font color="green">Sohrab Shah</font>*
  
  Lecture:
  
  Lab practical:
  
#### Links:

 * [PennCNV-Affy](http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/): In-depth guide into pre-processing of Affymetrix 6.0 microarrays for OncoSNP
 * [OncoSNP](https://sites.google.com/site/oncosnp/)
 * [Titan](http://compbio.bccrc.ca/software/titan/)
 * [SnpEff/SnpSift](http://snpeff.sourceforge.net/)


***

###  Module 6: Somatic Mutations <a id="module_6"></a>

  *<font color="green">Sohrab Shah</font>*
  
  Lecture:
  
  Lab practical:
  
#### Links:
 * [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)
 * [MutationSeq](http://compbio.bccrc.ca/software/mutationseq/)
  
  
***

##  Day 4 <a id="day_4"></a>

###  Module 7: Gene Expression Profiling <a id="module_7"></a>

  *<font color="green">Paul Boutros</font>*
  
  Lecture:
  
  Lab practical:
  
#### Links:
 * [R Statistical Package](http://www.r-project.org/) 
 * [R Studio](http://rstudio.org/)
  
  
***

###  Module 8: Variants to Networks <a id="module_8"></a>

#### Part 1: How to annotate variants and prioritize potentially relevant ones <a id="part_1"></a>
  
  *<font color="green">Juri Reimand</font>*
  
  Lecture:
  
  Lab practical:
  
  
***

#### Part 2: From genes to pathways <a id="part_2"></a>
  
  *<font color="green">Juri Reimand</font>*
  
  Lecture:
  
  Lab practical:
  

***

##  Day 5 <a id="day_5"></a>

#### Part 3: Network Analysis using Reactome <a id="part_3"></a>

*<font color="green">Lincoln Stein and Robin Haw</font>*

  Lecture:
  
  Lab practical:
  
#### Papers:

Clustering Algorithms: [Newman Clustering](http://www.pnas.org/content/103/23/8577.abstract) and [Hotnet](http://www.ncbi.nlm.nih.gov/pubmed/22174262)

Reactome Website: [NAR paper](http://www.ncbi.nlm.nih.gov/pubmed/26656494)

[Nature Methods and Perspectives Paper](http://www.ncbi.nlm.nih.gov/pubmed/26125594)

[Supplementary Materials](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4717906/bin/NIHMS750562-supplement-Supplementary_materials.pdf)


#### Links:

 * [GO](http://www.geneontology.org)
 * [KEGG](http://www.genome.jp/kegg)
 * [Biocarta](http://www.biocarta.com)
 * [Reactome](http://reactome.org/) Curated human pathways
 * [NCI/PID](http://pid.nci.nih.gov/)
 * [Pathway Commons](http://www.pathwaycommons.org/pc/) Aggregates pathways from multiple sources
 * [iRefWeb/iRefIndex](http://wodaklab.org/iRefWeb/) Protein interactions
 * [>300 more](http://www.pathguide.org/)

#### Tools for finding/converting gene identifiers and gene attributes

 * [Ensembl/BioMart](http://www.ensembl.org/index.html) 
 * [The Synergizer](http://llama.mshri.on.ca/synergizer/translate/) 

#### Cytoscape
 
 * [Cytoscape ](http://www.cytoscape.org/)
 * [Open Tutorials for Cytoscape](http://opentutorials.cgl.ucsf.edu/index.php/Portal:Cytoscape)
 
Useful plugins:

  * VistaClara - makes it easy to visualize gene expression data on networks
  * Agilent Literature Search - extracts interactions from PubMed abstracts
  * clusterMaker - provides multiple ways to cluster gene expression and networks
  * BiNGO - provides over-representation analysis using Gene Ontology in Cytoscape - you can select genes in your network or provide a list of genes and see the enrichment results visually mapped to the Gene Ontology
  * commandTool, coreCommands - used to control Cytoscape by a series of commands. E.g. automate the process: open network, layout network, save network as PDF. These plugins require Cytoscape 2.7
  * jActiveModules - requires gene expression data over multiple samples (>3). Finds regions of a network where genes are active (e.g. differentially expressed) across multiple samples.
  * [EnrichmentMap](http://baderlab.org/Software/EnrichmentMap)
  * [ReactomeFI](http://wiki.reactome.org/index.php/Reactome_FI_Cytoscape_Plugin)
  * [Many more](http://chianti.ucsd.edu/cyto_web/plugins/index.php)

  
  
***

###  Module 9: Integration of Clinical Data <a id="module_9"></a>

*<font color="green">Anna Lapuk</font>*

  Lecture:
  
  Lab practical:
  
#### Papers:

[Cox Regression Survival Paper](http://www.ncbi.nlm.nih.gov/pubmed/12613350)

[Genomics and Transcriptional Abberations Paper](http://www.ncbi.nlm.nih.gov/pubmed/17157792)

[Supplementary Data](http://www.cell.com/cms/attachment/588951/4507288/mmc1.pdf)

  
  
  ***
