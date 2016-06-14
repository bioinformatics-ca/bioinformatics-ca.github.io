---
layout: post2
permalink: /informatics_for_rna_seq_analysis_2016/
title: Informatics for RNA-Seq Analysis 2016 Student Page
header1: Informatics for RNA-Seq Analysis 2016
header2: Workshop pages for students
image: CBW_RNA_seq_icon.jpg
---

<ul id="navmenu">
  <li><a id="back_to_top">Contents</a>
     <ul class="sub1">
     <li><a href="#course_schedule">Course Schedule</a></li>
     <li><a href="#q_a_forum">Workshop Q/A Forum</a></li>
     <li><a href="#laptop_setup">Laptop Setup Instructions</a></li>
     <li><a href="#pre_readings">Pre-Workshop Readings</a></li>
     <li><a href="#pre_tutorials">Pre-Workshop Tutorials</a></li>
     <li><a href="#amazon_cloud">Amazon Cloud</a></li>
     <li><a href="#r_review">R Review Session</a></li>
      <li><a href="#day1">Day 1</a>
         <ul class="sub2">  
           <li><a href="#welcome">Welcome</a></li>
           <li><a href="#module_1">Module 1</a></li>
           <li><a href="#module_2">Module 2</a></li>
           <li><a href="#module_3">Module 3</a></li>
           <li><a href="#assignment1">Integrated Assignment</a></li>
        </ul>
      </li>
       <li><a href="#day_2">Day 2</a>
          <ul class="sub2">
             <li><a href="#module_4">Module 4</a></li>
             <li><a href="#module_5">Module 5</a></li>
           </ul>
       </li>
    </ul>
  </li>
</ul>  

<br>

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/rna-seq/RNA-Seq_2016_Schedule_v1.pdf">Schedule for June 16 to June 17, 2016</a>

[&uarr;](#back_to_top)

###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/RNASeq2016">here</a>!

[&uarr;](#back_to_top)

###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/rna-seq/laptop_setup_instructions.pdf">here</a>.

[&uarr;](#back_to_top)
  
##### Difference Between **R** and **RStudio**
<br>
**RStudio** doesn't know where libraries are installed, when they are not installed through the **RStudio** package manager. To tell **RStudio** the location, you can define the path in a startup file. Create a file called .Renviron . Inside there:

```r
R_LIBS=<R Library Path of other installed packages>
```

That was the problem when students installed things in **RStudio** at the command line using the **R** command <code>install.package()</code>.

... or you could use the package manger to install libraries.

[&uarr;](#back_to_top)

##### Syntax highlighting
<br>
... of scripts in the R editor does not seem to work under Windows. If you want highlighted syntax, use RStudio instead.

[&uarr;](#back_to_top)

###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop read these:
  
  [Integrative Genomics Viewer (IGV): high-performance genomics data visualization and exploration](http://www.ncbi.nlm.nih.gov/pubmed/22517427)
  
  [Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks](http://www.ncbi.nlm.nih.gov/pubmed/22383036)
  
  [ENCODE RNA-Seq Standards](https://genome.ucsc.edu/ENCODE/protocols/dataStandards/ENCODE_RNAseq_Standards_V1.0.pdf)
  
  [Methods to study splicing from high-throughput RNA sequencing data](http://www.ncbi.nlm.nih.gov/pubmed/24549677)
  
  [Differential analysis of gene regulation at transcript resolution with RNA-seq](http://www.ncbi.nlm.nih.gov/pubmed/23222703)
  
  [A comprehensive assessment of RNA-seq accuracy, reproducibility and information content by the Sequencing Quality Control Consortium](http://www.ncbi.nlm.nih.gov/pubmed/25150838)
  
  [Recurrent chimeric RNAs enriched in human prostate cancer identified by deep sequencing](http://www.ncbi.nlm.nih.gov/pubmed/21571633)
  
[&uarr;](#back_to_top)

###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [R Tutorial](http://www.cyclismo.org/tutorial/R/) up to and including 5. Basic Plots
* The [R command cheat sheet](../../resources/R_Short-refcard.pdf)

2) **UNIX Preparation tutorials**: 

* [UNIX Bootcamp](https://github.com/griffithlab/rnaseq_tutorial/wiki/Unix-Bootcamp)
* Tutorials #1-3 on [UNIX Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)
* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 

3) **IGV Tutorial**: Review how to use IGV Genome Browser if you have not used this tool before.

* The [IGV Tutorial](http://bioinformatics-ca.github.io/bioinformatics_for_cancer_genomics_IGV_lab_2016/)

4) [Sequencing Terminology](http://www.ncbi.nlm.nih.gov/projects/genome/glossary.shtml)

[&uarr;](#back_to_top)

### Logging into the Amazon Cloud <a id="amazon_cloud"></a>

Instructions can be found [here](http://bioinformatics-ca.github.io/logging_into_the_Amazon_cloud/).

* We have set up 32 instances on the Amazon cloud - one for each student. In order to log in to your instance, you will need a security certificate. If you plan on using Linux or Mac OS X, please download [this certificate](http://cbwmain.dyndns.info/private/CBWCG.pem). Otherwise if you plan on using Windows (with Putty and Winscp), please download [this certificate](http://cbwmain.dyndns.info/private/CBWCG.ppk).

[&uarr;](#back_to_top)

***

### R Review Session <a id="r_review"></a>

*<font color="827e9c">Fouad Yousif</font>*

[Lecture](../../resources/RReview_slides.pdf)

[Scripts](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/resources/R_Review_Session_Code.ipynb)

[&uarr;](#back_to_top)

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="827e9c">Ann Meyer</font>* 
<br>

[&uarr;](#back_to_top)

###  Module 1: Introduction to Cloud Computing <a id="module_1"></a>

  *<font color="827e9c">Obi Griffith</font>*
  
  Lecture:
  
  Lab practical:

[&uarr;](#back_to_top)

###  Module 2: Introduction to RNA Sequencing and Analysis <a id="module_2"></a>

  *<font color="827e9c">Malachi Griffith</font>*
  
  Lecture:
  
  Lab practical:

[&uarr;](#back_to_top)

###  Module 3: RNA-Seq Alignment and Visualization <a id="module_3"></a>

  *<font color="827e9c">Fouad Yousif</font>*
  
  Lecture:
  
  Lab practical

[&uarr;](#back_to_top)  
  
### Integrated Assignment <a id="assignment"></a>

*<font color="827e9c">Fouad Yousif</font>*

Paper: [Recurrent chimeric RNAs enriched in human prostate cancer identified by deep sequencing](http://www.ncbi.nlm.nih.gov/pubmed/21571633)

[Assignment Questions](https://github.com/griffithlab/rnaseq_tutorial/wiki/Integrated-Assignment)

Assignment Answers:

[&uarr;](#back_to_top)

##  Day 2 <a id="day_2"></a>

###  Module 4: Expression and Differential Expression <a id="module_4"></a>

  *<font color="827e9c">Obi Griffith</font>*
  
  Lecture:
  
  Lab practical:

[&uarr;](#back_to_top)

###  Module 5: Isoform Discovery and Alternative Expression <a id="module_5"></a>

  *<font color="827e9c">Malachi Griffith</font>*
  
  Lecture:
  
  Lab Practical:

[&uarr;](#back_to_top)

***

### Keeping Up-to-date with RNA-Seq Analysis Developments

For additional resources, tutorials, future directions, and more please refer to the [RNA-seq wiki](http://www.rnaseq.wiki/)

[&uarr;](#back_to_top)

***

