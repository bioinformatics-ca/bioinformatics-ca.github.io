---
layout: post2
permalink: /exploratory_analysis_of_biological_data_2016/
title: Exploratory Analysis of Biological Data 2016 Student Page
header1: Exploratory Analysis of Biological Data 2016
header2: Workshop pages for students
image: CBW_R_icon.jpg
---

<ul id="navmenu">
  <li><a id="back_to_top">Contents</a>
     <ul class="sub1">
     <li><a href="#course_schedule">Course Schedule</a></li>
     <li><a href="#q_a_forum">Workshop Q/A Forum</a></li>
     <li><a href="#setup">Laptop Setup Instructions</a></li>
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
             <li><a href="#module_5">Module 5</a></li>
             <li><a href="#assignment2">Integrated Assignment</a></li>
           </ul>
       </li>
    </ul>
  </li>
</ul>  

<br>

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/Stats/Stats_2016_Schedule_v1.pdf">Schedule for June 7 to June 8, 2016</a>

[&uarr;](#back_to_top)

###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://noteapp.com/EDAinR">here</a>!
  
### Workshop Survey

Your feedback is important to us.  Please complete our [workshop survey](https://www.surveymonkey.com/r/2C5BTPB).
  
### Class Photo

![Class photo](https://bioinformatics-ca.github.io/2016_workshops/Stats/CBW-June-7.jpeg)

[&uarr;](#back_to_top)

### Setup Instructions <a id="setup"></a>


Let's begin - start here by installing the first "project" we will work with, as soon as you come in:

1. Place a **pink** PostIt on your laptop screen so we know you are working on these instructions.
2. Make sure you have a "workshop directory" (perhaps called "training") in which you install all the RStudio "projects" we will be working with.
3. You must have **R** installed. If you don't, [**install R now**](http://cran.utstat.utoronto.ca/).
4. You must have **RStudio** installed. If you don't, [**install RStudio now**](https://www.rstudio.com/products/rstudio/download/).
5. open **RStudio**
6. Select **File &rarr; NewProject...**
7. Click on **Version Control**
8. Click on **Git**
9. Enter `https://github.com/hyginn/R_EDA-Introduction` as the **Repository URL**.
10. Click on **Browse...** to find your **training** directory...
11. (The *project-directory name* should autofill to `R_EDA-Introduction`)
12. Click **Create Project**; the project files should be downloaded and the *console* should prompt you to type `init()` to begin.
13. Type `init()` into the *console* pane.
14. Place a **green** PostIt if this has worked for you; place a **pink** PostIt if you run into issues or have questions.

---
<!--
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
* The [R command cheat sheet](../../resources/R_Short-refcard.pdf)

[&uarr;](#back_to_top)

***
-->

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="#827e9c">Ann Meyer</font>* 
<br>

[&uarr;](#back_to_top)


###  Module 1: Exploratory Data Analysis <a id="module_1"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  [Lecture](https://bioinformatics.ca/eda-module-1-2016)
  
---

**Before we move on: [let's list our progress!](https://docs.google.com/document/d/1l6vE8ImKPuE-s8Xltpx9NpSl6YaE4t2pY8qicH05Ihs/edit?usp=sharing)
**

---

  
<!--
  Resources:
  
  * [Single-cell RNA-Seq defines cell type](http://www.ncbi.nlm.nih.gov/pubmed/24531970)
  * [Beyond Bar and Line-Graphs](http://www.ncbi.nlm.nih.gov/pubmed/25901488)
-->  
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
  
[&uarr;](#back_to_top)

***


###  Module 2: Regression Analysis <a id="module_2"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  [Lecture](https://bioinformatics.ca/eda-module-2-2016)
  
---

**In RStudio, make a new project from `https://github.com/hyginn/R_EDA-Regression` as the Repository URL**

**Before we move on: [let's list our progress!](https://docs.google.com/document/d/1l6vE8ImKPuE-s8Xltpx9NpSl6YaE4t2pY8qicH05Ihs/edit?usp=sharing)
**

---
  
  Links:
  
* [Maximal Information Coefficient](http://www.ncbi.nlm.nih.gov/pubmed/22174245)
* [Homepage for data exploration with the MIC measure](http://www.exploredata.net/) 
* [**CRAN**: package MINERVA](http://cran.r-project.org/web/packages/minerva/)  (**R** wrapper for a fast *mine* implementation)

[&uarr;](#back_to_top)

***

###  Module 3: Dimension Reduction <a id="module_3"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  [Lecture](https://bioinformatics.ca/eda-module-3-2016)
  
---

**In RStudio, make a new project from `https://github.com/hyginn/R_EDA-DimensionReduction` as the Repository URL**

**Before we move on: [let's list our progress!](https://docs.google.com/document/d/1l6vE8ImKPuE-s8Xltpx9NpSl6YaE4t2pY8qicH05Ihs/edit?usp=sharing)
**

---

[&uarr;](#back_to_top)

***

### Integrated Assignment Part 1 <a id="assignment1"></a>

*<font color="#827e9c">Lauren Erdman and Ben Brew</font>*

[Lecture](https://bioinformatics.ca/eda-integrated-assignment-2016)

[Assignment](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/Stats/Stats2016_IntegratedAssignment.docx)

[Assignment Questions](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/Stats/Stats2016_IntegratedAssignment_Answers.R)

[Data Set](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/Stats/ccleCgc.rda)

[Plot.R](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/Stats/plot.R)

[&uarr;](#back_to_top)

***

##  Day 2 <a id="day_2"></a>

###  Module 4: Clustering Analysis <a id="module_4"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  [Lecture](https://bioinformatics.ca/eda-module-4-2016)
  
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

[&uarr;](#back_to_top)

###  Module 5: Hypothesis Testing for EDA <a id="module_5"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  [Lecture](https://bioinformatics.ca/eda-module-5-2016)
  
  Scripts:
  
  Links:
  
  * [NGS Differential Transcriptional Analysis](http://www.ncbi.nlm.nih.gov/pubmed/25894390)
  * [Erroneous Analysis of Significance](http://www.ncbi.nlm.nih.gov/pubmed/21878926)

[&uarr;](#back_to_top)

###  Integrated Assignment Part 2 <a id="assignment2"></a>

  *<font color="#827e9c">TBA</font>*
  
  Assignment Part 2:
  
  Questions in **R** Part 2:
  
  Answer Key Part 2:
  
  [&uarr;](#back_to_top)
  
  ***
