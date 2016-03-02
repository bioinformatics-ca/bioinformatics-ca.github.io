---
layout: post2
permalink: /informatics_for_rna_seq_analysis_2016/
title: Informatics for RNA-Seq Analysis 2016 Student Page
header1: Informatics for RNA-Seq Analysis 2016
header2: Workshop pages for students
image: CBW_RNA_seq_icon.jpg
---

#### Contents
[Course Schedule](#course_schedule)

[Workshop Q/A Forum](#q_a_forum)

[Laptop Setup Instructions](#laptop_setup)

[Pre-Workshop Readings](#pre_readings)

[Pre-Workshop Tutorials](#pre_tutorials)

[Logging into the Amazon Cloud](#amazon_cloud)

...[Logging in with ssh (Mac/Linux)](#ssh_login)

...[Logging in with Putty (Windows)](#putty_login)

...[File System Layout](#file_system_layout)

[R Review Session](#r_review)

**[Day 1](#day_1)**


  ...[Welcome](#welcome)
  
  ...[Module 1: Introduction to Cloud Computing](#module_1)
  
  ...[Module 2: Introduction to RNA Sequencing and Analysis](#module_2)
  
  ...[Module 3: RNA-Seq Alignment and Visualization](#module_3)
  
  ...[Integrated Assignment](#assignment)
  
  
**[Day 2](#day_2)**


  ...[Module 4: Expression and Differential Expression](#module_4)
  
  ...[Module 5: Isoform Discovery and Alternative Expression](#module_5)
  

***

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/rna-seq/RNA-Seq_2016_Schedule_v1.pdf">Schedule for June 16 to June 17, 2016</a>


###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/RNASeq2016">here</a>!


###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/rna-seq/laptop_setup_instructions.pdf">here</a>.
  
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

###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop read these:
  
  [Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration](http://www.ncbi.nlm.nih.gov/pubmed/22517427)
  
  [Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.ncbi.nlm.nih.gov/pubmed/22383036)
  
  [ENCODE RNA-Seq Standards](https://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf)
  
  [Methods to study splicing from high-throughput RNA sequencing data](http://www.ncbi.nlm.nih.gov/pubmed/24549677)
  
  [Differential analysis of gene regulation at transcript resolution with RNA-seq](http://www.ncbi.nlm.nih.gov/pubmed/23222703)
  
  [A comprehensive assessment of RNA-seq accuracy, reproducibility and information content by the Sequencing Quality Control Consortium](http://www.ncbi.nlm.nih.gov/pubmed/25150838)
  
  [Recurrent chimeric RNAs enriched in human prostate cancer identified by deep sequencing](http://www.ncbi.nlm.nih.gov/pubmed/21571633)
  
###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [R Tutorial](http://www.cyclismo.org/tutorial/R/) up to and including 5. Basic Plots
* The [R command cheat sheet](../../resources/R_Short-refcard.pdf)

2) **UNIX Preparation tutorials**: 

* [UNIX Bootcamp](https://github.com/griffithlab/rnaseq_tutorial/wiki/Unix-Bootcamp)
* Tutorials #1-3 on [UNIX Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)
* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 

3) **IGV Tutorial**: Review how to use IGV Genome Browser if you have not used this tool before.

* The [IGV Tutorial](../../resources/IGV_Tutorial.pdf)
* *Datasets too large for Github*

4) **Sequencing Terminology**:

* *Coming soon*

### Logging into the Amazon Cloud <a id="amazon_cloud"></a>

* These instructions will **ONLY** be relevant in class, as the Cloud will not be accessible from home in advance of the class.
 
* We have set up 30 instances on the Amazon cloud - one for each student. In order to log in to your instance, you will need a security certificate. If you plan on using Linux or Mac OS X, please download this **link certificate here**. Otherwise if you plan on using Windows (with Putty and Winscp), please download this **link certificate here**. 

* On the cloud, we're going to use the default username: **ubuntu**

#### Logging in with ssh (Mac/Linux) <a id="ssh_login"></a>
<br>
##### Logging in

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



##### Logging in



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

### R Review Session <a id="r_review"></a>

*<font color="827e9c">Fouad Yousif</font>*

Lecture:

Scripts:


##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="827e9c">Ann Meyer</font>* 
<br>

###  Module 1: Introduction to Cloud Computing <a id="module_1"></a>

  *<font color="827e9c">Obi Griffith</font>*
  
  Lecture:
  
  Lab practical:


###  Module 2: Introduction to RNA Sequencing and Analysis <a id="module_2"></a>

  *<font color="827e9c">Malachi Griffith</font>*
  
  Lecture:
  
  Lab practical:


###  Module 3: RNA-Seq Alignment and Visualization <a id="module_3"></a>

  *<font color="827e9c">Fouad Yousif</font>*
  
  Lecture:
  
  Lab practical
  
  
### Integrated Assignment <a id="assignment"></a>

*<font color="827e9c">Fouad Yousif</font>*

Paper: [Recurrent chimeric RNAs enriched in human prostate cancer identified by deep sequencing](http://www.ncbi.nlm.nih.gov/pubmed/21571633)

[Assignment Questions](https://github.com/griffithlab/rnaseq_tutorial/wiki/Integrated-Assignment)

Assignment Answers:


##  Day 2 <a id="day_2"></a>

###  Module 4: Expression and Differential Expression <a id="module_4"></a>

  *<font color="827e9c">Obi Griffith</font>*
  
  Lecture:
  
  Lab practical:


###  Module 5: Isoform Discovery and Alternative Expression <a id="module_5"></a>

  *<font color="827e9c">Malachi Griffith</font>*
  
  Lecture:
  
  Lab Practical:


***

### Keeping Up-to-date with RNA-Seq Analysis Developments

For additional resources, tutorials, future directions, and more please refer to the [RNA-seq wiki](http://www.rnaseq.wiki/)

***

