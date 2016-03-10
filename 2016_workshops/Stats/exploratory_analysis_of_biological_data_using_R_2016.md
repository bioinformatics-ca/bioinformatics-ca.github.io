---
layout: post2
permalink: /exploratory_analysis_of_biological_data_2016/
title: Exploratory Analysis of Biological Data 2016 Student Page
header1: Exploratory Analysis of Biological Data 2016
header2: Workshop pages for students
image: CBW_R_icon.jpg
---

<ul id="navmenu">
  <li><a href="contents">Contents</a>
     <ul class="sub1">
     <li><a href="#course_schedule">Course Schedule</a></li>
     <li><a href="#q_a_forum">Workshop Q/A Forum</a></li>
     <li><a href="#laptop_setup">Laptop Setup Instructions</a></li>
     <li><a href="#pre_readings">Pre-Workshop Tutorials</a></li>
      <li><a href="#day1">Day 1</a>
         <ul class="sub2">  
           <li><a href="#welcome">Welcome</a></li>
           <li><a href="#data_sets">Data Sets for Modules 1 - 5</a></li>
           <li><a href="#module_1">Module 1</a></li>
           <li><a href="#module_2">Module 2</a></li>
           <li><a href="#module_3">Module 3</a></li>
           <li><a href="#assignment1">Integrated Assignment</a></li>
        </ul>
      </li>
       <li><a href="#day_2">Day 2</a>
          <ul class="sub2">
             <li><a href="#module_4">Module 4</a></li>
             <li><a href="#module_5">Nodule 5</a></li>
             <li><a href="#assignment2">Integrated Assignment</a></li>
           </ul>
       </li>
    </ul>
  </li>
</ul>  

<br>

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/Stats/Stats_2016_Schedule_v1.pdf">Schedule for June 7 to June 8, 2016</a>


###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/EDA2016">here</a>!


###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/population/laptop_setup_instructions.pdf">here</a>.
  
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
... of scripts in the **R** editor does not seem to work under Windows. If you want highlighted syntax, use **RStudio** instead.


###  Pre-workshop Tutorials <a id="pre_readings"></a>

  You need to be familiar with the material covered in the Introduction to **R** tutorial, below. The tutorial should be very accessible even if you have never used R before.
  
* The [Introduction to **R** tutorial](http://steipe.biochemistry.utoronto.ca/abc/index.php/R_tutorial) 
* The R command cheat sheet **put cheat sheet here**


***

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="green">Ann Meyer</font>* 
<br>


### Data Sets for Modules 1 - 5 <a id="data_sets"></a>

**link datasets here**

###  Module 1: Exploratory Data Analysis <a id="module_1"></a>

  *<font color="green">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Resources:
  
  * [Single-cell RNA-Seq defines cell type](http://www.ncbi.nlm.nih.gov/pubmed/24531970)
  * [Beyond Bar and Line-Graphs](http://www.ncbi.nlm.nih.gov/pubmed/25901488)
  
  Helpful Links:
  
* [The **R** help mailing list](https://stat.ethz.ch/mailman/listinfo/r-help)
* [**Rseek**: the specialized search engine for **R** topics](http://rseek.org/)
* [**R** questions on stackoverflow](http://stackoverflow.com/questions/tagged/r)
* [The Comprehensive **R** Archive Network **CRAN**](http://cran.r-project.org/)
* [The **CRAN** task-view collection](http://cran.r-project.org/web/views/)
* [**Bioconductor** task views](http://www.bioconductor.org/packages/release/BiocViews.html)
* [Using Projects with R Studio](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects)
* [Software Carpentry](http://software-carpentry.org/)
* [Best Practices for Scientific Computing](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745)
* [Version control in R Studio](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)
  
***

###  Module 2: Regression Analysis <a id="module_2"></a>

  *<font color="green">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Links:
  
* [Maximal Information Coefficient](http://www.ncbi.nlm.nih.gov/pubmed/22174245)
* [Homepage for data exploration with the MIC measure](http://www.exploredata.net/) 
* [**CRAN**: package MINERVA](http://cran.r-project.org/web/packages/minerva/)  (**R** wrapper for a fast *mine* implementation)


***

###  Module 3: Dimension Reduction <a id="module_3"></a>

  *<font color="green">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:


***

### Integrated Assignment Part 1 <a id="assignment1"></a>

*<font color="green">TBA</font>*

Part 1:

Part 2:

Data Set:

For your reference:

***

##  Day 2 <a id="day_2"></a>

###  Module 4: Clustering Analysis <a id="module_4"></a>

  *<font color="green">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Links:
  
* [Comparison of Clustering Methods](http://www.ncbi.nlm.nih.gov/pubmed/19240124)
* [**R**-"task view": Cluster Analysis](http://cran.r-project.org/web/views/Cluster.html)  (and Finite Mixture Models)
  
  Dataset:
  
  If you load using **Gset.RData** do:
  
```r
  load("gset.RData")
```
on the command line.  (Check that 'gset' is actually lower case in the folder.  You might need a capital letter at the start.)

If you load using **Platf.RData** do:

```r
  load("platf.RData")
```

R object file: **GSE26922.rds** 

Read with:

```r
  gset <- readRDS("GSE26922.rds")
  # do not run the following line:
  gset <- gset [ [ idx ] ]
```

###  Module 5: Hypothesis Testing for EDA <a id="module_5"></a>

  *<font color="green">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Links:
  
  * [NGS Differential Transcriptional Analysis](http://www.ncbi.nlm.nih.gov/pubmed/25894390)
  * [Erroneous Analysis of Significance](http://www.ncbi.nlm.nih.gov/pubmed/21878926)


###  Integrated Assignment Part 2 <a id="assignment2"></a>

  *<font color="green">TBA</font>*
  
  Assignment Part 2:
  
  Questions in **R** Part 2:
  
  Answer Key Part 2:
  
  ***
