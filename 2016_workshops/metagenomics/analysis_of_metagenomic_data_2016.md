---
layout: post2
permalink: /analysis_of_metagenomic_data_2016/
title: Analysis of Metagenomic Data 2016 Student Page
header1: Analysis of Metagenomic Data 2016
header2: Workshop pages for students
image: CBW_Metagenome_icon.jpg
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
  
  ...[Module 1: Introduction to Metagenomics and Computing in the Cloud](#module_1)
  
  ...[Module 2: Marker Gene-based Analysis of Taxonomic Composition](#module_2)
  
  ...[Module 3: Introduction to PICRUSt](#module_3)

  ...[Integrated Assignment Part 1](#assignment1)
  
  
**[Day 2](#day_2)**


  ...[Module 4: Metagenomic Taxonomic Composition](#module_4)
  
  ...[Module 5: Metagenomic Functional Composition](#module_5)
  
  ...[Integrated Assignment Part 2](#assignment2)
  
  
**[Day 3](#day_3)**

  
  ...[Module 6: Metatranscriptomics](#module_6)
  
  ...[Module 7: Biomarker Selection](#module_7)
  

***

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/metagenomics/Metagenomics_2016_Schedule_v1.pdf">Schedule for June 22 to June 24, 2016</a>


###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/Metagenomics2016">here</a>!


###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/ht-seq/laptop_setup_instructions.pdf">here</a>.


###  Pre-workshop Tutorials <a id="pre_tutorials"></a>

1) **UNIX Preparation tutorials**: Please complete tutorials #1-3 on UNIX at http://www.ee.surrey.ac.uk/Teaching/Unix/

* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 

2) **Cytoscape Preparation tutorials**: Complete the [introductory tutorial to Cytoscape](http://opentutorials.cgl.ucsf.edu/index.php/Portal:Cytoscape3)

* Introduction to Cytoscape3 - User Interface

* Introduction to Cytoscape3 - Welcome Screen

* Filtering and Editing in Cytoscape 3


###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop, read these.
  
  * [Bioinformatics for the Human Microbiome Project](http://www.ncbi.nlm.nih.gov/pubmed/23209389)
  
  * [Microbiome science needs a healthy dose of scepticism](http://www.ncbi.nlm.nih.gov/pubmed/25143098)
  
  * [Methylotrophic methanogenic Thermoplasmata implicated in reduced methane emissions from bovine rumen](http://www.ncbi.nlm.nih.gov/pubmed/23385573)
  
  
### Logging into the Amazon Cloud <a id="amazon_cloud"></a>

* These instructions will **ONLY** be relevant in class, as the Cloud will not be accessible from home in advance of the class.

* We have set up 30 instances on the Amazon cloud - one for each student. In order to log in to your instance, you will need a security certificate. If you plan on using Linux or Mac OS X, please download this **link certificate here**. Otherwise if you plan on using Windows (with Putty and Winscp), please download this **link certificate here**. 

* On the cloud, we're going to use the default username: **ubuntu**

#### Logging in with ssh (Mac/Linux) <a id="ssh_login"></a> <br>


##### Logging in <br>


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

#### Logging in with Putty (Windows) <a id="putty_login"></a><br>

##### Logging in<br>

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

###  Module 1: Introduction to Metagenomics and Computing in the Cloud <a id="module_1"></a>

  *<font color="green">William Hsiao</font>*
  
  Lecture:
  
  Lab practical:


***

###  Module 2: Marker Gene-based Analysis of Taxonomic Composition <a id="module_2"></a>

  *<font color="green">William Hsiao</font>*
  
  Lecture:
  
  Lab practical:
  
***

###  Module 3: Introduction to PICRUSt <a id="module_3"></a>

  *<font color="green">Morgan Langille</font>*
  
  Lecture:
  
  Lab practical:


***

### Integrated Assignment Part 1<a id="assignment1"></a>

*<font color="green">Thea Van Rossum and Mike Peabody</font>*

Integrated Assignment: 


***

##  Day 2 <a id="day_2"></a>


###  Module 4: Metagenomic Taxonomic Composition <a id="module_4"></a>

  *<font color="green">Morgan Langille</font>*
  
  Lecture:
  
  Lab practical:
  
***

###  Module 5: Metagenomic Functional Composition <a id="module_5"></a>

  *<font color="green">Morgan Langille</font>*
  
  Lecture:
  
  Lab practical:
  
***

### Integrated Assignment Part 2<a id="assignment2"></a>

*<font color="green">Thea Van Rossum and Mike Peabody</font>*

Integrated Assignment: 


***

##  Day 3 <a id="day_3"></a>

###  Module 6: Metatranscriptomics  <a id="module_6"></a>

  *<font color="green">John Parkinson</font>*
  
  Lecture:
  
  Lab practical:

***

###  Module 7: Biomarker Selection  <a id="module_7"></a>

  *<font color="green">Fiona Brinkman</font>*
  
  Lecture:

***
