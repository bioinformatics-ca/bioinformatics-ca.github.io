---
layout: post2
permalink: /intror_module2-3_2016/
title: IntroR Module 2-3 2016
header1: IntroR Module 2-3 2016
header2: Workshop pages for students
image: CBW_introtoR-icon.jpg
---

 
 ================================================== 
 
 Canadian Bioinformatics Workshops Series           
 Toronto, May 20 2015                               
 Introduction to R                                  
                                                    
 Faculty: Boris Steipe <boris.steipe@utoronto.ca>   
                                                    
 Module 3: Data Analysis                            
                                                    
 ================================================== 




 
### "for" loop



convert columns to numeric

~~~r  
str(sup3)

str(sup3[,2])

sup3[ ,2] <- as.numeric(sup3[ ,2])

for (i in 1:10) {
	print(i*i)
}

for (i in 2:ncol(sup3)) {
	sup3[,i] <- as.numeric(sup3[ ,i])
}  
~~~





### Filtering data


 Task: find the top 30 most differentially enriched
 genes in the Mo cells. Hint: you will need to sort
 results ... but `sort()` is not the function you need,
 you need `order()`.

 Then find the same for Mf cells.

 Then find the union of the two sets. Plot the 
 differential enrichment of one against the other.






### Basic plots and slightly more advanced plots


 Task: Plot boxplots for the different cell-types,
       then plot the actual values of requested
       genes.


 Task: show the differential expression as a
       barplot.

 Barplots are bad. Improve according to 
 Weissgerber et al.'s ideas


 More plotting topics



 
### Integrating data


 Task:
  
 1 - For the top 10 Monocy.
 
 2 - translate their gene symbols to Entrez
     IDs using <http://biodbnet.abcc.ncifcrf.gov/>
     
 3 - see whether they are co-expressed (i.e.
     presumably coregulated) at <http://coxpresdb.jp/>



