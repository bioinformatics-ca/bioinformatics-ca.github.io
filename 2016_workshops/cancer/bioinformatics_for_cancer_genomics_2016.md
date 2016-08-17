---
layout: post2
permalink: /bioinformatics_for_cancer_genomics_2016/
title: Bioinformatics for Cancer Genomics 2016 Student Page
header1: Bioinformatics for Cancer Genomics 2016
header2: Workshop pages for students
image: CBW_cancerDNA_icon-16.jpg
---

<ul id="navmenu">
  <li><a id="back_to_top">Contents</a>
     <ul class="sub1">
     <li><a href="#course_schedule">Course Schedule</a></li>
     <li><a href="#q_a_forum">Workshop Q/A Forum</a></li>
     <li><a href="#laptop_setup">Laptop Setup Instructions</a></li>
     <li><a href="#pre_tutorials">Pre-Workshop Tutorials</a></li>
     <li><a href="#pre_readings">Pre-Workshop Readings</a></li>
     <li><a href="#amazon_cloud">Amazon Cloud</a></li>
      <li><a href="#day1">Day 1</a>
         <ul class="sub2">  
           <li><a href="#welcome">Welcome</a></li>
           <li><a href="#module_1">Module 1</a></li>
           <li><a href="#module_2.1">Module 2.1</a></li>
           <li><a href="#module_2.2">Module 2.2</a></li>
           <li><a href="#r_review">R Review Session</a></li>
        </ul>
      </li>
       <li><a href="#day_2">Day 2</a>
          <ul class="sub2">
             <li><a href="#module_3">Module 3</a></li>
             <li><a href="#module_4">Module 4</a></li>
           </ul>
       </li>
       <li><a href="#day_2">Day 3</a>
          <ul class="sub2">
             <li><a href="#module_5">Module 5</a></li>
             <li><a href="#module_6">Module 6</a></li>
           </ul>
       </li>
       <li><a href="#day_2">Day 4</a>
          <ul class="sub2">
             <li><a href="#module_7">Module 7</a></li>
             <li><a href="#module_8">Module 8</a>
             	<ul class="sub3">
             		<li><a href="#part_1">Part 1</a></li>
             		<li><a href="#part_2">Part 2</a></li>
           		</ul>
             </li>
           </ul>
       </li>
       <li><a href="#day_2">Day 5</a>
          <ul class="sub2">
             <li><a href="#module_8">Module 8</a>
             	<ul class="sub3">
             		<li><a href="#part_3">Part 3</a></li>
           		</ul>
             </li>
             <li><a href="#module_9">Module 9</a></li>
           </ul>
       </li>
    </ul>
  </li>
</ul>  

<br>

###  Course Schedule  <a id="course_schedule"></a>

  <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/BiCG_2016_Schedule_final.pdf">Schedule for May 30 to June 3, 2016</a>

[&uarr;](#back_to_top)

###  Workshop Q/A Forum <a id="q_a_forum"></a>

  Post your workshop questions <a href="https://noteapp.com/CancerGenomics2016">here</a>!

[&uarr;](#back_to_top)

### Workshop Survey

We appreciate your feedback on your experience at the workshop.  Please complete our [survey](https://www.surveymonkey.com/r/MWHBLCS) at the end of the workshop.

###  Laptop Setup Instructions <a id="laptop_setup"></a>

  Instructions to setup your laptop can be found <a href="http://bioinformatics-ca.github.io/2016_workshops/cancer/laptop_setup_instructions.pdf">here</a>.

[&uarr;](#back_to_top)
  
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

[&uarr;](#back_to_top)

###  Pre-Workshop Tutorials <a id="pre_tutorials"></a>

1) **R Preparation tutorials**: You are expected to have completed the following tutorials in **R** beforehand. The tutorial should be very accessible even if you have never used **R** before.

* The [CBW R tutorial](http://bioinformatics-ca.github.io/CBW_R_Tutorial/) or [R Tutorial](http://www.cyclismo.org/tutorial/R/) 
* The [R command cheat sheet](../../resources/R_Short-refcard.pdf)
* [R Plotting Reference](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/resources/Plotting.Reference.ipynb)


2) **Cytoscape 3.x Preparation tutorials**: Complete the [introductory tutorial to Cytoscape 3.x](http://opentutorials.cgl.ucsf.edu/index.php/Portal:Cytoscape3): 

* Introduction to Cytoscape3 - User Interface
* Introduction to Cytoscape3 - Welcome Screen
* Introduction to Cytoscape 3.1 - Networks, Data, Styles, Layouts and App Manager


3) **UNIX Preparation tutorials**: Please complete tutorials #1-3 on [UNIX Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)

* [Unix Cheat sheet](http://www.rain.org/~mkummel/unix.html) 

[&uarr;](#back_to_top)

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

[&uarr;](#back_to_top)  
  
### Logging into the Amazon Cloud <a id="amazon_cloud"></a>

Instructions can be found [here](http://bioinformatics-ca.github.io/logging_into_the_Amazon_cloud/).
 
* We have set up 30 instances on the Amazon cloud - one for each student. In order to log in to your instance, you will need a security certificate. If you plan on using Linux or Mac OS X, please download [this certificate](http://cbwmain.dyndns.info/private/CBWCG.pem). Otherwise if you plan on using Windows (with Putty and Winscp), please download [this certificate](http://cbwmain.dyndns.info/private/CBWCG.ppk).

[&uarr;](#back_to_top)


### Class Photo

![Class Picture](http://bioinformatics-ca.github.io/images/cbw-may-30-bicg.jpeg)
[Link to Download Class Photo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/cbw-may-30-bicg.jpeg)

***

##  Day 1 <a id="day_1"></a>

###  Welcome <a id="welcome"></a>

  *<font color="#827e9c">Ann Meyer</font>* 
<br>

[&uarr;](#back_to_top)

***

###  Module 1: Introduction to Cancer Genomics <a id="module_1"></a>

  *<font color="#827e9c">Trevor Pugh</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-1-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=9mKkQOf1Qxs)

[&uarr;](#back_to_top)

***

###  Module 2.1: Databases and Visualization Tools <a id="module_2.1"></a>

  *<font color="#827e9c">Michelle Brazas and Florence Cavalli</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-2-part-1-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=tKqw-G7AofE)
  
  [Lab practical for ICGC](https://bioinformatics.ca/bicg-module-2-lab-2016)
  
  [Lab practical for IGV](http://bioinformatics-ca.github.io/bioinformatics_for_cancer_genomics_IGV_lab_2016/)
  
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
  
[&uarr;](#back_to_top)

***  

###  Module 2.2: Logging into the Cloud <a id="module_2.2"></a>

  *<font color="#827e9c">Francis Ouellette</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-2-part-2-2016)

[&uarr;](#back_to_top)

***

###  *Optional* **R** Review Session <a id="r_review"></a>

*<font color="#827e9c">Florence Cavalli</font>*

[Lecture](https://bioinformatics.ca/bicg-module-r-review-2016-lecture)

[R commands](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BiCG_2016_Rreview_Code.R)

#### Links:

 * [R Statistical Package](http://www.r-project.org/) 
 * [R Studio](http://rstudio.org/) 

[&uarr;](#back_to_top)

***

##  Day 2 <a id="day_2"></a>

###  Module 3: Mapping and Genome Rearrangement <a id="module_3"></a>

  *<font color="#827e9c">Jared Simpson</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-3-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=1i8mDfsMh4s)
  
  Lab practicals: [Part 1 - Mapping](mapping) and [Part 2 - Rearrangements](rearrangement)
  
#### Links:

 * [What does my SAM flag mean?](https://broadinstitute.github.io/picard/explain-flags.html)
 * [Tools for Mapping High-throughput Sequencing Data Paper (2012)](http://www.ncbi.nlm.nih.gov/pubmed/23060614) 
 * [SAM/BAM file specifications](http://samtools.sourceforge.net/SAM1.pdf) 
 * [samtools](http://samtools.sourceforge.net/) 
 * [bwa](http://bio-bwa.sourceforge.net/) 
 * [lumpy-sv](https://github.com/arq5x/lumpy-sv) 
  
[&uarr;](#back_to_top)

***

###  Module 4: Gene Fusion Discovery <a id="module_4"></a>

  *<font color="#827e9c">Andrew McPherson</font>*
  
[Lecture and Lab](https://bioinformatics.ca/bicg-module-4-2016)
  
Lab practical:

 * [Setup instructions](/2016_workshops/cancer/gene_fusions/install.html)
 * [Part 1 - Prediction](/2016_workshops/cancer/gene_fusions/run.html)
 * [Part 2 - Exploration](/2016_workshops/cancer/gene_fusions/exploration.html)
 * [Part 3 - Visualization](/2016_workshops/cancer/gene_fusions/visualization.html)

#### Papers and Background Material:

 * [A survey of best practices for RNA-seq data analysis](http://dx.doi.org/10.1186/s13059-016-0881-8)
 * [The impact of translocations and gene fusions on cancer causation](http://dx.doi.org/10.1038/nrc2091)
 * [The emerging complexity of gene fusions in cancer](http://dx.doi.org/10.1038/nrc3947)
 * [The landscape and therapeutic relevance of cancer-associated transcript fusions](http://dx.doi.org/10.1038/onc.2014.406)
 * [Fusion genes and their discovery using high throughput sequencing](http://dx.doi.org/10.1016/j.canlet.2013.01.011)

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

[&uarr;](#back_to_top)

***

##  Day 3 <a id="day_3"></a>

###  Module 5: Copy Number Alterations <a id="module_5"></a>

*<font color="#827e9c">Sohrab Shah and Fong Chun Chan</font>*
  
[Lecture](https://bioinformatics.ca/bicg-module-5-2016)
  
[Lab practical](https://bioinformatics.ca/bicg-module-5-lab-2016) 

  * [Lab Module](/2016_workshops/cancer/cnas/cna_lab.html)
    + This is the instructions for the lab practical.

  * [Data Analysis Package](/2016_workshops/cancer/cnas/cna_data_analysis_package.tar.gz)
    + Contains the various files and Rmarkdown file that will be used to do further exploration and analysis on copy number alterations.
    + This is package is already on the server. You can also download this to your own computer and perform the analyses locally.

  * [Software Installation](/2016_workshops/cancer/cnas/cna_installation.html)
    + This page contains information on how to install the different software used in the lab practical.

  * [Data Preparation](/2016_workshops/cancer/cnas/cna_prepdata.html)
    + This page contains information on how the data was prepared to be used for lab practical.

#### Data for Lab Practical

* [METABRIC Seg File](/2016_workshops/cancer/cnas/segs/METABRIC_DatasetI997.seg)
    + Seg file from the METABRIC project to be visualized in IGV.

#### Plots for Lab Practical

These plots are provided for convenience. They can be generated by following the lab practical.

  * Oncosnp
    + [HCC1395 OncoSNP Plot for Ploidy Configuration 1](/2016_workshops/cancer/cnas/plots/HCC1395.1.pdf)
    + [HCC1395 OncoSNP Plot for Ploidy Configuration 2](/2016_workshops/cancer/cnas/plots/HCC1395.2.pdf)

  * Titan
    + [HCC1395 TITAN Plot](/2016_workshops/cancer/cnas/plots/HCC1395_exome_tumour.all.pdf)

#### Links:

 * [PennCNV-Affy](http://penncnv.openbioinformatics.org/en/latest/user-guide/affy/): In-depth guide into pre-processing of Affymetrix 6.0 microarrays for OncoSNP
 * [OncoSNP](https://sites.google.com/site/oncosnp/)
 * [Titan](http://compbio.bccrc.ca/software/titan/)
 * [SnpEff/SnpSift](http://snpeff.sourceforge.net/)

[&uarr;](#back_to_top)

***

###  Module 6: Somatic Mutations <a id="module_6"></a>

  *<font color="#827e9c">Sohrab Shah and Fong Chun Chan</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-6-2016)
  
   <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=mKP39iTn1JQ)

  [Lab practical](https://bioinformatics.ca/bicg-module-6-lab-2016) 

  * [Lab Module](/2016_workshops/cancer/snvs/snv_lab.html)
    + This is the instructions for the lab practical.

  * [Data Analysis Package](/2016_workshops/cancer/snvs/snv_data_analysis_package.tar.gz)
    + Contains the various files and Rmarkdown file that will be used to do further exploration and analysis on somatic mutations data.
    + This is package is already on the server. You can also download this to your own computer and perform the analyses locally.

  * [Data Preparation](/2016_workshops/cancer/snvs/snv_prepdata.html)
    + This page contains information on how the data was prepared to be used for lab practical.

  * [Pre-processing Bams](/2016_workshops/cancer/snvs/snv_prepro_bam.html)
    + This page contains information on how to pre-process bam (e.g. filtering) for downstream analyses.
  
#### Links:

 * [Strelka](https://sites.google.com/site/strelkasomaticvariantcaller/)
 * [MutationSeq](http://compbio.bccrc.ca/software/mutationseq/)
  
[&uarr;](#back_to_top)
  
***

##  Day 4 <a id="day_4"></a>

###  Module 7: Gene Expression Profiling <a id="module_7"></a>

  *<font color="#827e9c">Fouad Yousif</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-7-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=m-QoR6QES9A)
  
  [Lab practical](https://bioinformatics.ca/bicg-module-7-questions-2016) [with answers](https://bioinformatics-ca.github.io/2016_workshops/cancer/assignment_questions_with_answers.pdf)
  
#### Links:
 * [R Statistical Package](http://www.r-project.org/) 
 * [R Studio](http://rstudio.org/)
  
[&uarr;](#back_to_top)
  
***

###  Module 8: Variants to Networks <a id="module_8"></a>

#### Part 1: How to annotate variants and prioritize potentially relevant ones <a id="part_1"></a>
  
  *<font color="#827e9c">Robin Haw</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-8-part-1-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=Ch7OvLxjDl8)
  
  [Lab practical](https://bioinformatics-ca.github.io/BiCG_Module8_Annovar_lab/)
  
  [Data Set Input - VCF](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part1_Reimand/Passed.somatic.snvs.vcf)
  
  [Data Set Output - Annovar text table](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part1_Reimand/passed.somatic.snvs.vcf.annovar.out.txt.hg19_multianno.txt)
  
  [Lab Practical extra info](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part1_Reimand/BiCG_2016_Module8_Part1_VariantAnn_AnnovarInstall2016.txt)
  
#### Links
  
  [Annovar](http://annovar.openbioinformatics.org/en/latest/)
  
[&uarr;](#back_to_top)
  
***

#### Part 2: From genes to pathways <a id="part_2"></a>
  
  *<font color="#827e9c">Juri Reimand</font>*
  
  [Lecture](https://bioinformatics.ca/bicg-module-8-part-2-2016)
  
  [Lab practical protocol](https://bioinformatics-ca.github.io/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/BICG_2016_Module8_Part2_Protocol.pdf)
  
  Data Sets Gene Lists:
  
  [Data Set Genelist GBM](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/Genelist_GBM.txt)
  
  [Data Set Genelist KIRC](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/Genelist_KIRC.txt)
  
  Data Sets Enrichment Results (g:Profiler) from Gene Lists:
  
  [gProfiler Results GBM](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/gprofiler_results_GBM.txt)
  
  [gprofiler Results KIRC](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/gprofiler_results_KIRC.txt)
  
  [gProfiler hsapiens](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/hsapiens.NAME.gmt)
  
  Data Sets Enrichment Map (Cytoscape) from Enrichment Results:
  
  [EM cys](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/LAB_cytoscape_session.cys)
  
  Enrichmentmap
  
  [Enrichment Map](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/BICG_2016_Module8-Part2_Reimand/enrichmentmap-2.0.1.jar)

[&uarr;](#back_to_top)  

***

##  Day 5 <a id="day_5"></a>

#### Part 3: Network Analysis using Reactome <a id="part_3"></a>

*<font color="#827e9c">Robin Haw</font>*

  [Lecture](https://bioinformatics.ca/bicg-module-8-part-3-2016)
  
  <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=ae1YkHDDG7E)
  
  [Lab practical](https://bioinformatics.ca/bicg-module-8-part-3-lab-2016) and [Answers](https://bioinformatics-ca.github.io/2016_workshops/cancer/BiCG_2016_Module8_LabAnswers-v2.pdf)
  
  Data Sets:
  
  [OVCA_TCGA_Clinical.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/OVCA_TCGA_Clinical.txt)
  
  [OVCA_TCGA_GeneList.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/OVCA_TCGA_GeneList.txt)
  
  [OVCA_TCGA_MAF.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/OVCA_TCGA_MAF.txt)
  
#### Papers:

[Integrated genomic analyses of ovarian carcinoma](http://www.nature.com/nature/journal/v474/n7353/full/nature10166.html)

Clustering Algorithms: [Newman Clustering](http://www.pnas.org/content/103/23/8577.abstract) and [Hotnet](http://www.ncbi.nlm.nih.gov/pubmed/22174262)

Reactome Website: [NAR paper](http://www.ncbi.nlm.nih.gov/pubmed/26656494); [Website guide](https://bioinformatics-ca.github.io/2016_workshops/cancer/ReactomeWebSiteGuide_for_resources.pdf)

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

[&uarr;](#back_to_top)  
  
***

###  Module 9: Integration of Clinical Data <a id="module_9"></a>

*<font color="#827e9c">Anna Goldenberg and Lauren Erdman</font>*

  [Lecture](https://bioinformatics.ca/bicg-module-9-2016)
  
   <img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/images/Nova-Video-300px.png?raw=true" width="42"> [Recorded Lecture](https://www.youtube.com/watch?v=ysuwIjj5KR0)
  
  [Updated Lab](https://bioinformatics.ca/bicg-module-9-lab-2016)
  
  [RData](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/Module9/OICR-Survival-Workshop-Data-revised-6-3-2016.RData)
  
  [R Commands](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/Module9/OICR%20Data%20Integration%20and%20Survival%20Workshop%20Script-6-3-2016.R)
  
#### Tools:

[Predict tool](http://www.predict.nhs.uk/predict.html)
  
#### Papers:

[Gene expression profiling reveals molecularly and clinically distinct subtypes of glioblastoma multiforme](http://www.pnas.org/content/102/16/5814.full)

[Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1](http://www.ncbi.nlm.nih.gov/pubmed/20129251)

[Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis](http://www.ncbi.nlm.nih.gov/pubmed/19759197)

[Similarity network fusion for aggregating data types on a genomic scale](http://www.nature.com/nmeth/journal/v11/n3/abs/nmeth.2810.html)

 [&uarr;](#back_to_top) 
  
  ***

## Data for the Workshop ##

### Tool Installation ###

Instructions for installing the tools used in the workshops can be found [here](http://bioinformatics-ca.github.io/install_tools_2016/).

### Data Sets ###
- HCC1395 data: [CEL](http://www.hpc4health.ca/cbw/2016/CG_data/HCC1395/cel.tar.gz) [exome](http://www.hpc4health.ca/cbw/2016/CG_data/HCC1395/exome.tar.gz) [rnaseq](http://www.hpc4health.ca/cbw/2016/CG_data/HCC1395/rnaseq.tar.gz)
- [Module 3 data](http://www.hpc4health.ca/cbw/2016/CG_data/Module3.tar.gz)
- Module 4 data: [bams](http://www.hpc4health.ca/cbw/2016/CG_data/Module4/bams.tar.gz), [cbw_tutorial](http://www.hpc4health.ca/cbw/2016/CG_data/Module4/cbw_tutorial.tar.gz), [refdata](http://www.hpc4health.ca/cbw/2016/CG_data/Module4/refdata.tar.gz), [sampledata](http://www.hpc4health.ca/cbw/2016/CG_data/Module4/sampledata.tar.gz)
- Module 5 data: [data](http://www.hpc4health.ca/cbw/2016/CG_data/Module5.tar.gz) [ref_data](http://www.hpc4health.ca/cbw/2016/CG_data/ref_data.tar.gz)
- [Module 6 data](http://www.hpc4health.ca/cbw/2016/CG_data/Module6.tar.gz)
- [Module 7 data](http://www.hpc4health.ca/cbw/2016/CG_data/Module7.tar.gz)
 
### Results from Instructor's Instance on Amazon ###
- [Module3 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module3_result.tar.gz)
- [Module4 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module4/analysis.tar.gz)
- [Module5 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module5_result.tar.gz)
- [Module6 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module6_result.tar.gz)
- [Module7 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module7_result.tar.gz)
- [Module8 result](http://www.hpc4health.ca/cbw/2016/CG_data/Module8_result.tar.gz)
