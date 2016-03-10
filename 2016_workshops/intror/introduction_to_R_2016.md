---
layout: post2
permalink: /introduction_to_r_2016/
title: Introduction to R 2016 Student Page
header1: Introduction to R 2016
header2: Workshop pages for students
image: CBW_introtoR-icon.jpg
---

<ul id="navmenu">
  <li><a id="back_to_top">Contents</a>
     <ul class="sub1">
     <li><a href="#course_schedule">Course Schedule</a></li>
     <li><a href="#q_a_forum">Workshop Q/A Forum</a></li>
     <li><a href="#laptop_setup">Laptop Setup Instructions</a></li>
     <li><a href="#helpful_materials">Helpful Materials</a></li>
     <li><a href="#pre_readings">Pre-Workshop Readings</a></li>
      <li><a href="#day1">Day 1</a>
         <ul class="sub2">  
           <li><a href="#welcome">Welcome</a></li>
           <li><a href="#module_1">Module 1</a></li>
           <li><a href="#module_2">Module 2</a></li>
           <li><a href="#module_3">Module 3</a></li>
        </ul>
      </li>
    </ul>
  </li>
</ul>  

<br>

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/intror/IntroR_2016_Schedule_v1.pdf">Schedule for June 6, 2016</a>

[&uarr;](#back_to_top)

###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="http://todaysmeet.com/IntroR2016">here</a>!

[&uarr;](#back_to_top)

###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/intror/laptop_setup_instructions.pdf">here</a>.

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
... of scripts in the **R** editor does not seem to work under Windows. If you want highlighted syntax, use **RStudio** instead.

[&uarr;](#back_to_top)

###  Pre-workshop Readings <a id="pre_readings"></a>

  Before coming to the workshop, read these.

[&uarr;](#back_to_top)  
  
### Helpful Materials <a id="helpful_materials"></a>

* The [Introduction to **R** tutorial](http://steipe.biochemistry.utoronto.ca/abc/index.php/R_tutorial) 
* The R command cheat sheet **put cheat sheet here**

[&uarr;](#back_to_top)

***

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="#827e9c">Ann Meyer</font>* 
<br>

[&uarr;](#back_to_top)

###  Module 1: The **R** Environment <a id="module_1"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  * **R** script template
  * First steps
  
  Data:
  
  Practical:
  
* [Using Projects with R Studio](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects)
* [Software Carpentry](http://software-carpentry.org/)
* [Best Practices for Scientific Computing](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745)
* [Version control in R Studio](https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN)
  
  Resources:
  
  * [Single-cell RNA-Seq defines cell type](http://www.ncbi.nlm.nih.gov/pubmed/24531970)
  
  Helpful Links:
  
* [The **R** help mailing list](https://stat.ethz.ch/mailman/listinfo/r-help)
* [**Rseek**: the specialized search engine for **R** topics](http://rseek.org/)
* [**R** questions on stackoverflow](http://stackoverflow.com/questions/tagged/r)
* [The Comprehensive **R** Archive Network **CRAN**](http://cran.r-project.org/)
* [The **CRAN** task-view collection](http://cran.r-project.org/web/views/)
* [**Bioconductor** task views](http://www.bioconductor.org/packages/release/BiocViews.html)
  
[&uarr;](#back_to_top)
  
***

###  Module 2: Programming Basics <a id="module_2"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Data:
  
  Practical:

[&uarr;](#back_to_top)

***

###  Module 3: Using **R** for Data Analysis <a id="module_3"></a>

  *<font color="#827e9c">Boris Steipe</font>*
  
  Lecture:
  
  Scripts:
  
  Resources:
  
  * [Beyond Bar and Line-Graphs](http://www.ncbi.nlm.nih.gov/pubmed/25901488)
  
  Practical:
  
  * [Biological identifier conversion]( http://biodbnet.abcc.ncifcrf.gov/)
  * [Gene co-expression database](http://coxpresdb.jp/) 

[&uarr;](#back_to_top)

***
