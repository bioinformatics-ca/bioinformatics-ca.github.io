---
layout: post2
permalink: /PNAOD_module4_lab_2016/
title: Pathway and Network Analysis of -Omics Data 2016 Module 4 Lab
header1: Pathway and Network Analysis of -Omics Data 2016
header2: Module 4 Lab
image: CBW_pathway_icon.jpg
---

**This work is licensed under a [Creative Commons Attribution-ShareAlike 3.0 Unported License](http://creativecommons.org/licenses/by-sa/3.0/deed.en_US). This means that you are able to copy, share and modify the work, as long as the result is distributed under the same license.**

# *De Novo* Subnetwork Clustering Analysis: Reactome

By Robin Haw

## Aim

This exercise will provide you with an opportunity to perform and network analysis using the Reactome Functional Interaction (FI) and the ReactomeFIViz app.

## Goal 

Analyze somatic mutation data to identify biology that contributes to ovarian cancer.

### Example 1: Network-based analysis of OvCa somatic mutation data 

*	Open up Cytoscape. 

*	Go to Apps>Reactome FI and Select “Gene Set/Mutational Analysis”.

*	Choose “2015 (Latest)” Version. 

*	Upload/Browse [OVCA_TCGA_MAF.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/OVCA_TCGA_MAF.txt) file. 

*	Select “NCI MAF” (Mutation Annotation File) and Choose sample cutoff value of 4. 

*	**Do not** select “Fetch FI annotations”. 

*	Click OK.

1.	Describe the size and composition of the OvCa network?

2.	What are the most frequently mutated genes?

3.	Describe the TP53-PEG3 interaction, and the source information to support this interaction?

4.	Describe the data sources for the RB1-SMARCA4 FI?

5.	After clustering, how many modules are there? 

6.	How many pathway gene sets are there in Module 0 when the FDR Filter is set to 0.005 and Module Size Filter to 10?

*	Hint: Analyze Module Functions>Pathway Enrichment. Select appropriate filters at each step.

7.	What are the most significant pathway gene sets in Module 0, 2, 3 and 5? 

8.	Do the GO Biological Process annotations correlate with the significant pathway annotations for Module 0? 

*	Hint: Analyze Module Functions>GO Biological Process. Select appropriate filters at each step.

9.	What are the most significant GO Cell Component gene sets in Module 3? [Optional]

*	Hint: Analyze Module Functions>GO Cell Component. Select appropriate filters at each step.

10.	Are any of the modules annotated with the NCI Disease term: “Stage_IV_Breast_Cancer” [malignant cancer]?

*	Hint: Load Cancer Gene Index>Neoplasm>Neoplasm_by_Site>Breast Neoplasm>…….

11.	How many modules are statistically significant in the CoxPH analysis? 

*	Hint: Analyze Module Functions>Survival Analysis>Upload/Browse [OVCA_TCGA_Clinical.txt](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/cancer/OVCA_TCGA_Clinical.txt). Click OK.

12.	What does the Kaplan-Meyer plot show for the most clinically significant modules?

*	Hint: Click the most statistically significant module link [blue line] from the CoxPH results panel. Click OK. Click #_plot.pdf to display Kaplan-Meyer plot. Repeat this for the other significant module links. KM plot: samples having genes mutated in a module (red line), and samples having no genes mutated in the module (green line).

13.	Taking into what you have learned about module 2, what is your hypothesis?
