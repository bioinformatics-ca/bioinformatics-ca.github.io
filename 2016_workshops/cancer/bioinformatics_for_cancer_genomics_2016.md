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
  
##### Difference Between **R and **RStudio

RStudio doesn't know where libraries are installed, when they are not installed through the RStudio package manager. To tell RStudio the location, you can define the path in a startup file. Create a file called .Renviron . Inside there:

'''R
R_LIBS=<R Library Path of other installed packages>
'''

That was the problem when students installed things in **RStudio at the command line using the **R command <code>install.package()</code>.

... or you could use the package manger to install libraries.

##### Syntax highlighting

... of scripts in the R editor does not seem to work under Windows. If you want highlighted syntax, use RStudio instead.


###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

Do these before coming to the workshop.


###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop, read these.
  
  
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
