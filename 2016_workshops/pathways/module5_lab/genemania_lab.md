---
layout: post2
permalink: /PNAOD_module5_lab_genemania_2016/
title: Pathway and Network Analysis of -Omics Data 2016 Module 5 Lab
header1: Pathway and Network Analysis of -Omics Data 2016
header2: Module 5 Lab
image: CBW_pathway_icon.jpg
---

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# Module 5 Practical Lab: GeneMANIA (web version)

By Quaid Morris and Veronique Voisin 


## Goal of this practical lab 

Create GeneMANIA networks starting from a single gene to predict its function or starting from a gene list. Explore and understand the main output features of GeneMANIA such as the network composition or the enriched functions. 

This practical consists of 3 exercises. You can choose to do these exercises using the questions as your only guide - or see the following pages for the step-by-step checklist to finding the answers. 

Before starting the exercises,download the files:

*	[30_prostate_cancer_genes.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/30_prostate_cancer_genes.txt)

*	[Mixed_gene_list.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/mixed_gene_list.txt)

*	[CYP11B_pearson_correlation_prostate.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/CYB11B_pearson_correlation_prostate.txt)

Optional exercise: try to redo these exercises using the Cytoscape GeneMANIA application.

NOTE: Network layouts are flexible and can be rearranged. What you see when you perform these exercises may not be identical to what you see in the tutorial, or what you have seen other times that you have performed the exercises. Exact layouts and predictions can also be affected by updates to the networks database that GeneMANIA uses. However it is expected that the network weights and predicted genes will be similar to those shown here. 



## EXERCISE 1

1)	Imagine that you are interested in exploring the function of the human GRN gene: GRN returned as the strongest hit from your omics experiment but not many information about this gene is available in functional databases. Use GeneMANIA to identify its predicted function as well as potential interaction partners. 

**Skills**: GeneMANIA Single Gene search; Navigating Search Results; Exploring  available Genes features;  Rerun a new analysis using a single gene or multiple genes query from the network.

### STEPS

1. Go to GeneMANIA’s homepage at <http://www.genemania.org/>  

2. In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true).  

3. Enter the following gene: GRN	

4. Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true) and wait for the results.	

5. When your search results load, examine the network. Genes you searched with are indicated with stripes, related genes added by GeneMANIA are represented in black, and colored links represent the interactions that connect the nodes (genes). **Tip**: zoom in and zoom out using trackpad or mouse scrolling up and down using the mouse or trackpad or zoom in to the maximum using the specific icon ![arrow](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/arrow.png?raw=true).
 	
6. Clicking on a node gives information about the name of gene, the possibility to add, remove this gene or search with this gene only. Click on the GRN node and explore the displayed information. 	

7. Locate the Functions summary tab (bottom left icon ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)). What are the functions significantly associated with  this network?  GRN is the central node of this network: which function would you predict for GRN? How well did GeneMANIA perform (hints: use GeneCards (<http://www.genecards.org/>) , PubMed (<http://www.ncbi.nlm.nih.gov/pubmed/>))?	

8. Locate the gene with the strongest association with GRN. **Tip**: this gene is the largest node on the network. 	

9. Re-run an analysis by adding SORT1, HSPG2 to the search. Click on SORT1, click on Add ![add](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/add.png?raw=true). Do the same for HSPG2.  The genes are now added to the search box and a new network is automatically created.  Which functions are associated with this new network ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)? 	 

10. On the left side of the window are located icons that we haven’t yet explored. The first 3 buttons are activating different network layouts. Try the circular ![circular](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circledot.png?raw=true), the aligned ![aligned](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/twodown.png?raw=true), and the force_directed ![force](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/crossing.png?raw=true) layouts.  	

11. Choose your favorite layout and save the  network as an image  using the *Network image As Shown* option from the *save*  menu ![save](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/save.png?raw=true).  	

### EXERCISE 1 - STEPS 1-4


* Go to GeneMANIA’s homepage at <http://www.genemania.org/>.	

* In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true).	

* Enter the following gene: GRN.	

* Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true) and wait for the results.	

![GM1](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM1.png?raw=true) 

### EXERCISE 1 - STEP 5

* When your search results load, examine the network. Genes you searched with are indicated with stripes, related genes added by GeneMANIA are represented in black, and colored links represent the interactions that connect the nodes (genes). **Tip**: zoom in and zoom out using trackpad or mouse scrolling up and down using the mouse or trackpad or zoom in to the maximum using the specific icon ![arrow](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/arrow.png?raw=true). 	

![GM2](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM2.png?raw=true) 
 





### EXERCISE 1 - STEP 6


* Clicking on a node gives information about the name of gene, the possibility to add, remove this gene or search with this gene only. Click on the GRN node and explore the displayed information. 	

![GM3](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM3.png?raw=true)  

### Exercise 1 - STEP 7


* Locate the Functions summary tab (bottom left icon  ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)). What are the functions significantly associated with  this network ?  GRN is the central node of this network: which function would you predict for GRN?
How well did GeneMANIA perform (hints: use GeneCards (http://www.genecards.org/) , pubMed (http://www.ncbi.nlm.nih.gov/pubmed/) )?	

![GM4](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM4.png?raw=true)  

### Exercise 1 - STEP 8 


* Locate the gene with the strongest association with GRN. **Tip**: this gene is the largest node on the network. 


**Answer is SORT1**

### Execerise 1 - STEP 9


* Re-run an analysis by adding SORT1, HSPG2 to the search. Click on SORT1, click on Add ![add](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/add.png?raw=true)   . Do the same for HSPG2.  The genes are now added to the search box and a new network is automatically created.  Which functions are associated with this new network ( ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true))? 	 

![GM5](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM5.png?raw=true)

![GM6](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod6/GM6.png?raw=true)   


 













### Exercise 1 - STEP 10


* On the left side of the window are located icons that we haven’t yet explored. The first 3 buttons are activating different network layouts. Try the circular ![circular](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circledot.png?raw=true), the aligned ![aligned](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/twodown.png?raw=true) and the force_directed ![force](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/crossing.png?raw=true) layouts.  	

#### Circular layout

![GM7](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM7.png?raw=true) 
 
#### Aligned layout
 
![GM8](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM8.png?raw=true) 
 
#### Force directed layout
 
![GM9](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM9.png?raw=true) 

### Exercise 1 - STEP 11


* Choose your favorite layout and save the  network as an image  using the “Network image As Shown” option from the *save*  menu ![save](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/save.png?raw=true).  	

![GM10](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM10.png?raw=true)  


**Notes** about biological interpretation of the results: 

The top functions predicted by GeneMANIA for GRN were related to lysosome and vacuole.  A pubmed search could confirm these results: “We experimentally verified that granulin precursor (GRN) gene, whose mutations cause frontotemporal lobar degeneration, is involved in lysosome function.” (Transcriptional gene network inference from a massive dataset elucidates transcriptome organization and gene function. Belcastro et al. Nucleic Acids Res. 2011 Nov 1;39(20):8677-88. 2011. PMID:21785136)

**A paper describing the interaction between GRN and SORT1 and demonstrates how finding related genes could be relevant for elaborating therapy:**
 
Targeted manipulation of the sortilin–progranulin axis rescues progranulin haploinsufficiency. Lee et al. Hum Mol Genet. 2014 March 15; 23(6): 1467–1478. Published online 2013 October 26. doi: 10.1093/hmg/ddt534. PMCID:PMC3929086
“Progranulin (GRN) mutations causing haploinsufficiency are a major cause of frontotemporal lobar degeneration (FTLD-TDP). Recent discoveries demonstrating sortilin (SORT1) is a neuronal receptor for PGRN endocytosis and a determinant of plasma PGRN levels portend the development of enhancers targeting the SORT1–PGRN axis. We demonstrate the preclinical efficacy of several approaches through which impairing PGRN's interaction with SORT1 restores extracellular PGRN levels. “

![GM11](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM11.png?raw=true) 

---

## EXERCISE 2

To start this exercise, you need to download the [30_prostate_cancer_genes.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/30_prostate_cancer_genes.txt) file from the course wiki and save it on your computer. 

For this exercise, you are working with a list of 30 prostate cancer genes. This list can be downloaded after the workshop from the cBioPortal website (<http://www.cbioportal.org/>) under the download section, user defined list. The cBioPortal for Cancer Genomics stores genomic data from large scale, integrated cancer genomic data sets. During this exercise, you will explore the types of networks that have been used to create the GeneMANIA network from the prostate cancer gene list and you will see how changing input parameters can affect the results. The last step of the exercise consists of uploading a custom network which is a list of genes that are positively correlated with CYP11B1 in mRNA expression data of 94 prostate cancer samples (<http://www.cbioportal.org/>) .

**Skills**: GeneMANIA search using a gene list; Navigating Search Results; Exploring  Networks and advanced options; Uploading a custom network. 

### STEPS

1.	Go to GeneMANIA’s homepage at <http://www.genemania.org/>	

2.	In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true)  .	

3.	Copy and paste genes in the file [30_prostate_cancer_genes.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/30_prostate_cancer_genes.txt).  	

4.	Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true) and wait for the results.	

5.	When your search results load, examine the network. Genes you searched with are indicated with stripes, related genes added by GeneMANIA are represented in black, and colored links represent the interactions that connect the nodes (genes). Move nodes around by selecting them with a mouse to investigate how they are connected.
	
6.	Click any link (edge) connecting two nodes to highlight information about it. **Note**: Clicking on an edge between 2 nodes will display information about all interaction networks that connect  these 2 nodes. It indicates the reference (publication) for these interactions. The colors indicate the type of interaction (co-expression, shared protein domains, co-localization, physical interactions and predicted). 
	
7.	Locate the Networks summary tab (on the right ![lines](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/threelines.png?raw=true)) and look at what data has been used to create the network and predictions. **Note**that Co-expression (purple colored lines, weight over 30%) and Shared protein domains (lightgold colored  lines, weight over 30%) influence the results the most, but Co-localization (blue colored lines), Physical interactions (salmon colored lines) and Predicted (orange) data are also included. At the top of the Networks summary tab, try Expand “none”, then “top” and “all”  to get information about the sources of the different networks. 
	
8.	Highlight all connections corresponding to each network by clicking the name of each network category. Click on  “Shared protein domains” and see which genes are connected by predicted protein protein interaction. You can do the same for “Co-localization” , “Co-expression”  and “Physical interactions”. **Tips**:  these observations of the number of connections make easier to understand why co-expression and shared protein domains get the higher percent weight on this network: they are helping to connect more genes than physical interactions and predicted; A higher weight means that this network helped more to find related genes. 	

9.	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at what functions were  significantly enriched in this list of prostate genes.	

10.	“Shared protein domains” is an important part of the network. What would be the GeneMANIA results if we don’t include this source when we run GSEA? Click on ‘Show advanced option ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true)’ which is located at the right of the search box. Uncheck ‘Shared protein domains’ and click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)  . Explore the results.
	
11.	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at what functions were  significantly enriched with these new settings.	

12.	Upload a custom network to GeneMANIA: in ‘advanced options’, click on “Upload network…” and browse your computer to locate and select the file [CYP11B_pearson_correlation_prostate.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/CYB11B_pearson_correlation_prostate.txt). Wait about a minute for the network to be uploaded.  Explore the results and locate the genes linked by the custom network (**tip**: click on “Uploaded” in the Networks tab). 
	
13.	Try additional parameters of the ‘advanced options ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true)’ by selecting “Customise advanced options" and look at how the changes you made influenced the results. For example change ‘Network weighting’ method or ‘Max resultant genes: ’.  	



### Exercise 2 - STEPS 1 to 4


1.	Go to GeneMANIA’s homepage at <http://www.genemania.org/>	

2.	In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true).	

3.	Copy and paste genes in the file [30_prostate_cancer_genes.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/30_prostate_cancer_genes.txt).  	

4.	Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)    and wait for the results.	

![GM12](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM12.png?raw=true)  














### Exercise 2 - STEP 5


*	When your search results load, examine the network. Genes you searched with are indicated with stripes, related genes added by GeneMANIA are represented in black, and colored links represent the interactions that connect the nodes (genes). Move nodes around by selecting them with a mouse to investigate how they are connected.
	
![GM13](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM13.png?raw=true) 
 

### Exercise 2 - STEP 6.


*	Click any link (edge) connecting two nodes to highlight information about it. **Note**: Clicking on an edge between 2 nodes will display information about all interaction networks that connect  these 2 nodes. It indicates the reference (publication) for these interactions. The colors indicate the type of interaction (co-expression, shared protein domains, co-localization, physical interactions and predicted). 
	

![GM14](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM14.png?raw=true)  



### Exercise 2 - STEP 7


*	Locate the Networks summary tab (on the right ![lines](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/threelines.png?raw=true)) and look at what data has been used to create the network and predictions. **Note** that Co-expression (purple colored lines, weight over 30%) and Shared protein domains (lightgold colored  lines, weight over 30%) influence the results the most, but Co-localization (blue colored lines), Physical interactions (salmon colored lines) and Predicted (orange) data are also included. At the top of the Networks summary tab, try Expand “none”, then “top” and “all”  to get information about the sources of the different networks. 
	
![GM15](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM15.png?raw=true) 

![GM16](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM16.png?raw=true) 

 



 




### Exercise 2 - STEP 8

![GM17](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM17.png?raw=true) 

*	Highlight all connections corresponding to each network by clicking the name of each network category. Click on  “Shared protein domains” and see which genes are connected by predicted protein protein interaction. You can do the same for “Co-localization” , “Co-expression”  and “Physical interactions”. **Tips**:  these observations of the number of connections make easier to understand why co-expression and shared protein domains get the higher percent weight on this network: they are helping to connect more genes than physical interactions and predicted; A higher weight means that this network helped more to find related genes. 	

![GM18](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM18.png?raw=true) 
 

### Exercise 2 - STEP 9

*	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at what functions were  significantly enriched in this list of prostate genes.	

![GM19](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM19.png?raw=true) 
 
### Exercise 2 - STEP 10


*	“Shared protein domains” is an important part of the network. What would be the GeneMANIA results if we don’t include this source when we run GSEA? Click on show advanced option ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true) which is located at the right of the search box. Uncheck ‘Shared protein domains’ and click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)  . Explore the results.
	
![GM20](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM20.png?raw=true) 
 
![GM21](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM21.png?raw=true) 

 

### Exercise 2 - STEP 11


*	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at what functions were  significantly enriched with these new settings.	

![GM22](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM22.png?raw=true)  

### Exercise - STEP 12


*	Upload a custom network to GeneMANIA: in ‘advanced options ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true)’, click on “Upload network…” and browse your computer to locate and select the file [CYP11B_pearson_correlation_prostate.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/CYB11B_pearson_correlation_prostate.txt). Wait about a minute for the network to be uploaded.  Explore the results and locate the genes linked by the custom network (**tip**: click on “Uploaded” in the Networks tab). 
	
![GM23](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM23.png?raw=true) 

![GM24](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM24.png?raw=true) 

![GM25](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM25.png?raw=true) 
 

 



 

### Exercise 2 - STEP 13.


*	Try additional parameters of the ‘advanced options ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true) ’ by selecting “Customise advanced options" and look at how the changes you made influenced the results. For example change ‘Network weighting’ method or ‘Max resultant genes: ’.  	


![GM26](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM26.png?raw=true) 
 
**How to get the 30 prostate cancer gene list**:

![GM27](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM27.png?raw=true)   

---

## EXERCISE 3

To start this exercise, you need to download the [Mixed_gene_list.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/mixed_gene_list.txt) file and save it on your computer. 

For this exercise, you are working on a gene list created by combining 3 user defined gene lists available from the cBioportal (<http://www.cbioportal.org>). It contains genes implicated in the DNA damage response, the PI3K-AKT-mTOR signaling pathway and Folate transport. This list is representative of a gene list obtained from transcriptomics data. During this exercise, we will first characterize our gene list based on functions and then we  will add potential drug and microRNAs targeting genes in the network, and we will save the report.


**Skills**: GeneMANIA search using a gene list; Navigating Search Results; Exploring  Functions; Adding attributes; Create a report. 

### STEPS  

1.	Go to GeneMANIA’s homepage at <http://www.genemania.org/>.	

2.	In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true)  .	

3.	Copy and paste genes in the file [Mixed_gene_list.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/mixed_gene_list.txt). Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)    and wait for the results.	

4.	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at functions returned by GeneMANIA	

5.	In the functions summary tab,  check some functions to color genes included in these functions. To follow this tutorial, you can for example color the “response to insulin” , “DNA recombination” and“vitamin transport” functions. **Tip**: You need to scroll down to found DNA recombination and vitamin transport as these pathways have an FDR greater than the one for “response to insulin”. 	

6.	We are going next to add miRs and drug interaction networks. Click on ‘Show advanced option   ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true)’ which is located at the right of the search box. Check “Drug-interactions-2013” and “miRNA-target-predictions-2013” as “Attributes”. Check “Physical interactions” and “Co-expression” . Click on “Customise advanced options”. Set “Max resultant genes” to 20 and “Max resultant attributes” to 40. Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)    and wait for the results. Explore the network. **Tips**: the Drug-interactions and miRNA-target-predictions nodes are displayed in gray. The nodes connected to a drug are genes that are targeted by the drug and nodes connected to a  microRNA (miR) are genes predicted to be targeted by this miR.	

7.	Locate our favorite gene PDPK1 on the network, select it by moving the mouse cursor to its node and wait there for a second. It will highlight this gene and all its connections. **Tip**: make sure that the black box displaying the details for PDPK1 is closed.	

8.	Generate and save a report of your results by locating the save menu ![save](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/save.png?raw=true), and selecting “Report”. The PDF report provides a detailed description of your search and results.	

9.	Investigate the “history” function by clicking on the related icon ![redo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/redo.png?raw=true) located at the bottom of the window. A panel pops up showing the past networks generated by GeneMANIA. Clicking on one panel will relaunch the search for this network. 	



### Exercise 3 - STEPS 1 - 3

1.	Go to GeneMANIA’s homepage at <http://www.genemania.org/>.	

2.	In the search window, ensure that the model organism is set to *Homo sapiens* ![homo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Up.png?raw=true).	

3.	Copy and paste genes in the file [Mixed_gene_list.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/pathways/module5_lab/mixed_gene_list.txt). Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)    and wait for the results.	

![GM28](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM28.png?raw=true) 

![GM29](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM29.png?raw=true)  

 

### Exercise 3 - STEP 4

*	Locate the Functions summary tab (bottom left ![circle](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/circle.png?raw=true)) and look at functions returned by GeneMANIA	

![GM30](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM30.png?raw=true)  


### Exercise 3 - STEP 5

*	In the functions summary tab,  check some functions to color genes included in these functions. To follow this tutorial, you can for example color the “response to insulin” , “DNA recombination” and“vitamin transport” functions. **Tip**: You need to scroll down to found DNA recombination and vitamin transport as these pathways have an FDR greater than the one for “response to insulin”. 	

![GM31](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM31.png?raw=true)  


### Exercise 3 - STEPS 6 

*	We are going next to add miRs and drug interaction networks. Click on ‘Show advanced option ![options](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/dotdotdot.png?raw=true)   ’ which is located at the right of the search box. Check “Drug-interactions-2013” and “miRNA-target-predictions-2013” as “Attributes”. Check “Physical interactions” and “Co-expression” . Click on “Customise advanced options”. Set “Max resultant genes” to 20 and “Max resultant attributes” to 40. Click on the search icon ![search](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/Search.png?raw=true)    and wait for the results. Explore the network.

**Tips**: the Drug-interactions and miRNA-target-predictions nodes are displayed in gray. The nodes connected to a drug are genes that are targeted by the drug and nodes connected to a  microRNA (miR) are genes predicted to be targeted by this miR.	

![GM32](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM32.png?raw=true) 

![GM33](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM33.png?raw=true)  

 


### Exercise 3 - STEP 7

*	Locate our favorite gene PDPK1 on the network, select it by moving the mouse cursor to its node and wait there for a second. It will highlight this gene and all its connections. 
**Tip**: make sure that the black box displaying the details for PDPK1 is closed.	

![GM34](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM34.png?raw=true)  

### Exercise 3 - STEP 8

*	Generate and save a report of your results by locating the save menu ![save](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/save.png?raw=true), and selecting “Report”. The PDF report provides a detailed description of your search and results.	

![GM35](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM35.png?raw=true) 

![GM36](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM36.png?raw=true)  

 


### Exercise 3 - STEP 9

*	Investigate the “history” function by clicking on the related icon ![redo](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/redo.png?raw=true) located at the bottom of the window. A panel pops up showing the past networks generated by GeneMANIA. Clicking on one panel will relaunch the search for this network. 	

![GM37](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/pathways/img/mod5/GM37.png?raw=true)  




Definitions:

What are the networks: Definition of the types of interaction:

*	Shared domains: Protein domain data. Two gene products are linked if they have the same protein domain. These data are collected from domain databases, such as InterPro, SMART and Pfam.

*	Co-localization: Genes expressed in the same tissue, or proteins found in the same location. Two genes are linked if they are both expressed in the same tissue or if their gene products are both identified in the same cellular location.

*	Co-expression: Gene expression data. Two genes are linked if their expression levels are similar across conditions in a gene expression study. Most of these data are collected from the Gene Expression Omnibus (GEO); we only collect data associated with a publication.

*	Predicted: Predicted functional relationships between genes, often protein interactions. A major source of predicted data is mapping known functional relationships from another organism via orthology.

What is defined by evidence sources?:

*	Evidence sources are the information contained in the multiple databases that GeneMANIA uses to establish interaction between two genes.

Network:

*	Node : circle representing the genes

*	Edge: line that links two nodes and represent an interaction between two genes (multiple lines correspond to multiple sources

*	Node size: Mapped to gene score, i.e. the degree to which GeneMANIA predicts the genes are related

*	Thickness of edges: Strength/weight of interaction

Layout: The layout is different each time so the user can request the layout run multiple times until the user is satisfied with the result.

in Networks tab:

*	Percent weight (score): a higher weight means that this network helped more to find related genes.


in Functions tab:

*	FDR: False discovery rate (FDR) is greater than or equal to the probability that this is a false positive.

*	Coverage: (number of genes in the network with a given function) / (all genes in the genome with the function)

#### In advanced options: 

*	Network weighting? GeneMANIA can use a few different methods to weight networks when combining all networks to form the final composite network that results from a search. The default settings are usually appropriate, but you can choose a weighting method in the advanced option panel. (more details at <http://pages.genemania.org/help/>).

*	Related genes: are genes added by GeneMANIA in addition to the genes from the query. It helps to grow the network and then to predict function of the query gene(s).

*	The attributes represent  the differences sources of evidence that can be used to build the network.


**Notes**:

*	prostate cancer gene list is “AKR1C3 AR CYB5A CYP11A1 CYP11B1 CYP11B2 CYP17A1 CYP19A1 CYP21A2 HSD17B1 HSD17B10 HSD17B11 HSD17B12 HSD17B13 HSD17B14 HSD17B2 HSD17B3 HSD17B4 HSD17B6 HSD17B7 HSD17B8 HSD3B1 HSD3B2 HSD3B7 RDH5 SHBG SRD5A1 SRD5A3 STAR”.

*	mixed gene list is AKT1 AKT1S1 AKT2 ATM ATR BRCA1 BRCA2 CHEK1 CHEK2 FANCF FOLR1 FOLR2 FOLR3 FOXO1 FOXO3 MDC1 MLH1 MLST8 MSH2 MTOR PARP1 PDPK1 PIK3CA PIK3R1 PIK3R2 PTEN RAD51 RHEB RICTOR RPTOR SLC19A1 TSC1 TSC2

*	**Tip**: look at GeneMANIA help pages when you run an analysis on your own after the workshop: <http://pages.genemania.org/help/>.

