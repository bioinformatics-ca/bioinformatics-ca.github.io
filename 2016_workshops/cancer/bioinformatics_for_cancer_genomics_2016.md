---
layout: post
permalink: /bioinformatics_for_cancer_genomics_2016/
title: Bioinformatics for Cancer Genomics 2016 Student Page
header1: Bioinformatics for Cancer Genomics 2016
header2: Workshop pages for students
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
  
  ...[Module 2: Databases and Visualization Tools](#module_2)


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
  
  
  **[Day 4](#day_4)**
  
  
  ......[Part 3: How to annotate variants and prioritize potentially relevant ones ](#part_3)
  
  ...[Module 9: Integration of Clinical Data ](#module_9)
  
  

***

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/BiCG_2016_Schedule_v1.pdf">Schedule for May 30 to June 3, 2016</a>


###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/CancerGenomics2016">here</a>!


###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/laptop_setup_instructions.pdf">here</a>.
  
##### Difference Between **R** and **RStudio**


**RStudio** doesn't know where libraries are installed, when they are not installed through the **RStudio** package manager. To tell **RStudio** the location, you can define the path in a startup file. Create a file called .Renviron . Inside there:

```r
R_LIBS=<R Library Path of other installed packages>
```

That was the problem when students installed things in **RStudio** at the command line using the **R** command <code>install.package()</code>.

... or you could use the package manger to install libraries.

##### Syntax highlighting


... of scripts in the R editor does not seem to work under Windows. If you want highlighted syntax, use RStudio instead.


###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [CBW R tutorial](http://bioinformatics.ca/workshop_wiki/index.php/R_tutorial) or http://www.cyclismo.org/tutorial/R/ 
* The [R command cheat sheet](http://bioinformatics.ca/workshop_wiki/images/4/4c/Short-refcard.pdf)
* [PlottingReference.pdf](http://bioinformatics.ca/workshop_wiki/images/d/dc/PlottingReference.pdf)


2) **Cytoscape 3.x Preparation tutorials**: Complete the introductory tutorial to Cytoscape 3.x: http://opentutorials.cgl.ucsf.edu/index.php/Portal:Cytoscape3
* Introduction to Cytoscape3 - User Interface
* Introduction to Cytoscape3 - Welcome Screen
* Introduction to Cytoscape 3.1 - Networks, Data, Styles, Layouts and App Manager


3) **UNIX Preparation tutorials**: Please complete tutorials #1-3 on UNIX at http://www.ee.surrey.ac.uk/Teaching/Unix/
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

#### Logging in with ssh (Mac/Linux) <a id="ssh_login"></a>

#### Logging in with Putty (Windows) <a id="putty_login"></a>

#### File System Layout <a id="file_system_layout"></a>




***

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="green">Ann Meyer</font>* 
<br>

###  Module 1: Concepts in Population Genomics <a id="module_1"></a>

  *<font color="green">Philip Awadalla</font>*
  
  Lecture:
  
  Lab practical:


###  Module 2: Extensions in Population Genomics <a id="module_2"></a>

  *<font color="green">Philip Awadalla</font>*
  
  Lecture:


###  Module 3: Basic Concepts in QTLs <a id="module_3"></a>

  *<font color="green">Philip Awadalla and Stephen Montgomery</font>*
  
  Lecture:
  
  Lab practical:


##  Day 2 <a id="day_2"></a>

###  Module 4: Advanced Concepts in QTLs <a id="module_4"></a>

  *<font color="green">Philip Awadalla and Stephen Montgomery</font>*
  
  Lecture:
  
  Lab practical:


###  Module 5: Gene x Environment <a id="module_5"></a>

  *<font color="green">Philip Awadalla and Stephen Montgomery</font>*
  
  Lecture:


###  Module 6: Functional Annotation of Mutations and Haplotypes <a id="module_6"></a>

  *<font color="green">Philip Awadalla</font>*
  
  Lecture:
  
  Lab practical:
