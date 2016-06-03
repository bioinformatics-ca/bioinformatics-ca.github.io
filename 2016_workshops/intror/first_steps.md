---
layout: post2
permalink: /intror_first_steps_2016/
title: IntroR First Steps 2016
header1: IntroR First Steps 2016
header2: Workshop pages for students
image: CBW_introtoR-icon.jpg
---

 ==================================================  
 
 Canadian Bioinformatics Workshops Series           
 Toronto, June 6, 2016                                
 Introduction to R                                  
                                                    
 Faculty: Boris Steipe <boris.steipe@utoronto.ca>   
                                                    
 Module 1: The R Environment                        
                                                    
 ==================================================


 
### Environment and user interface


#### Set up your workspace


 Let's set the current "working directory" to your CBW
 workshop directory. You have a workshop directory, right?

 Where are you now?

~~~r
getwd()
~~~

 syntax to set the working directory...

~~~r
setwd("/path/of/your/directory")  #the argument is a string
setwd("/path/of/your/course\ directory")  #"escape" blank spaces
list.files()
~~~

 A note for Windows: the `\` (backslash) has a particular
 meaning in strings: it "escapes" the following character.
 Therefore something like `setwd("C:\My\R\files\")` will NOT
 work. You have to "escape the escape character" to turn it
 into a "literal" backslash: `setwd("C:\\My\\R\\files\\")` or
 R will translate for you: `setwd("C:/My/R/files/")`

 I usually write `setwd("whatever/path")` as the first command
 in my scripts. More on that when we start writing programs.

 But how do you know what the full path for the working
 directory is?

 There is a neat trick to get the path for a file (on the Mac).
 Drag and drop a file into a script. This will put the full path
 of the file into your script (... as the argument to a source()
 command). Does this work on Windows?
 Linux?

 Try this.
 But unfortunately this does not work in R-Studio. There, you can
 use the menu. Then copy the command into your script.




#### Using scripts.


 All your R work should always go into a script.


 Put the cursor in any line, or select a block of code
 then press
    <command><return> (Mac)
    <ctrl><r> (Win)
 to execute it in the console.

 Try:

~~~r
length(dir("~", pattern = ".txt"))
~~~

 In the console, use `<up-arrow>`, `<left-arrow>` etc. to
 retrieve and edit previous commands.

 You can use `source("filename")` to execute an entire
 script at once.



### Help, documentation other information


 Recapitulate ...

~~~r
?dir     #help("dir")
??dir    #same as help.search("dir")

apropos("^dir")

 # the "sos" package to look for functions everywhere
if (!require(sos, quietly=TRUE)) {
  install.packages("sos")
  library(sos)
}

ls("package:sos")          contents of "sos"
browseVignettes("sos")     documentation

# use "sos" to find more functions...
findFn("directory")    
~~~



### Getting data into R


#### Assigning variables


~~~r
x <- "1"
x
phi <- 2 * pi
phi
x <- pi
~~~

 Question: how many significant digits does R process?
 How do you print them all?

 Task: find out how to control the number of digits printed in a print() expression.

~~~r
sprintf() # gives the most control
sprintf("%50.49f", pi)  
~~~


#### Digression: Anatomy of a function ...


 Various functions exist to display the properties of R objects. Here
 is a function that combines them:

~~~r
typeInfo <- function(x) {
    print(x, digits=22)  
    cat("str:    ")                
    str(x)  
    cat("mode:   ", mode(x), "\n")
    cat("typeof: ", typeof(x), "\n")
    cat("class:  ", class(x), "\n")
     if there are attributes, print them too
    if (! is.null(attributes(x))) {
        cat("attributes:\n")
        print(attributes(x))
    }
}
~~~

 That's a useful utility to have. Let's take it apart.
 Now: where do we put it?

### Creating vectors


 Recapitulate:
 
~~~r  
v <- c(1, 1, 2, 3, 5, 8)
v
v <- c(v, 13, 25)
v

1:3
seq(-0.5, 0.5, by = 0.1)
rep("Ha", 3)

genes <- c("Spic", "Cebpb", "Lyz2", "Sfpi1", "Nfkbiz")  
~~~

 These are some of the genes that are markers for
 monocytes...

 Often our data can be copied and pasted: open the text file
 for [Fig_3-CharacteristicGenes](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/intror/Fig_3-CharacteristicGenes.txt)

 Task: how do we get this data into a vector?
       First think of a way to do this by hand.
       Then, design a function that would do this conveniently.


 Task:
 Consolidate working with text: convert a binomial scientific name into a 5-letter label.


#### Glueing vectors together to make matrices:


 Task: make a vector with cell-types, according to Fig. 3, use `rep()` expressions, don't type out every line.


 Task: Make a matrix from the gene names, so that each row contains the cell type and a characterisitc gene.


 Task: use `cbind()` to assemble the two vectors into one matrix.




## Lists


### Data frames
 

 A data frame is a matrix or a "set" of data. It is a list of vectors and/or factors of the same length, that are related "across", such that data in the same position come from the same experimental unit (subject, animal, etc).

 Importantly, the columns can have different type!

~~~r  
myDF <- data.frame(genes = c("Abc1", "Qrz", "Fubr31"),
                   expr = c(168059, 23490578, 34),
                   induced = c(TRUE, FALSE, FALSE))
myDF[-2,]
typeInfo(myDF)  
~~~

 What is it with the factors? ...


~~~r  
myDF <- data.frame(genes = c("Abc1", "Qrz", "Fubr31"),
                   expr = c(168059, 23490578, 34),
                   induced = c(TRUE, FALSE, FALSE),
                   stringsAsFactors = FALSE)
typeInfo(myDF)  
~~~


 Why do we need data frames if they do the much the same as a list?
 More efficient storage, and indexing! 
 R's read...() functions return data frames.

 ... which gets us back to our task of looking at expression values.


 The relevant data is in supplementary table 3 which is an Excel spreadsheet.

 Download it from the workshop Wiki:



 A word on Excel: it's a very good spreadsheet program, it is miserable and often wrong on statistics, and it makes horrible, horrible plots.

 To elaborate - see the two links below:
 
 <http://www.practicalstats.com/xlsstats/excelstats.html>
 
 <http://www.burns-stat.com/documents/tutorials/spreadsheet-addiction/>
 
 ... these are not merely cosmetic problems!

 Therefore: Ok to keep data in Excel spreadsheets if you must - but read it into R for any analysis!

 But be cautious: one of the problems of Excel is that it truncates numeric precision. 
 Protip: convert all cells to "text" before export.

 There are many other "read" functions.
 Refer to the R Data Import / Export manual
 http://cran.r-project.org/doc/manuals/R-data.html
 See:
 
~~~r
?read.table  #... includes read.csv and read.delim
?read.fwf    #... for "fixed width format"
?readLines   #... for reading in text-files line by line
~~~

 Excel spreadsheets should be converted to csv or tsv format. Alternatively the package xlsreadwrite is available via CRAN ... see
 <http://cran.r-project.org/web/packages/xlsReadWrite/> ... but I think this is unsound practice.

 Task:
 
 1 - load the data in [Table_S3.xls](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/raw/master/2016_workshops/intror/Table_S3.xls) into Excel, and save it as a .csv (comma separated values) file.
 
 2 - Examine the file (a plain text file) in a text-editor (such as R). 
 
 3 - Read the table into R, assigning it to a variable. I usually give the first input of data the variable name "rawDat" since it will usually be processed before it becomes meaningful for analysis.
 
 4 - Use head() to look at the beginning of the object.
 
 5 - Remove any unneeded header rows.
 
 6 - Give the columns names that reflect the cell type (cf.  Figure 2c), and the stimulus status.
 
 7 - Use typeInfo() to analyse the object you have created.

 Much output. For a heavy-duty function, we should rewrite type info to limit the output ...

 Now: what is it with the "factors".


### Digression: Factors...


 Many of R's dataframe methods convert all strings into factors by default. Factors are special types: they are nominally strings - (m, f) or (agree, neutral, disagree) or such. But underlyingly they are coded as integers that identify the "levels" of the factor.

 To illustrate.
 
~~~r
genders <- factor(c("m", "f", "f", "m", "f"))
genders
typeInfo(genders)
is.ordered(genders)
~~~

 We can define ordered factors by adding some details to
 our factor() command - e.g. tumor grades:

~~~r
sampleGrades <- factor(c("G1", "G3", "G3", "G1", "G2", "G1"),
                       levels = c("G1", "G2", "G3", "G4"),
                       ordered = TRUE)
sampleGrades    #Note that G4 is a level although it was not observed in the data
is.ordered(sampleGrades) 
~~~

 Factors are useful since they support a number of analysis methods such as ordering boxplots, or calculating 

 For more on factors, have a look at this factor tutorial by Jenny Bryan: 
 http://www.stat.ubc.ca/~jenny/STAT545A/block08_bossYourFactors.html
 and this discussion on their use:
 http://stackoverflow.com/questions/3445316/factors-in-r-more-than-an-annoyance
 

 But for our purposes, the default behavior of R, to treat all strings as factors is entirely unwanted and needs to be turned off. Always use the parameter stringsAsFactors = FALSE to achieve this. If you don't you are in for some bad surprises if e.g. there is a character "contaminant" in a numeric column.

~~~r
myDF <- data.frame(data = c("N/A", 1, 1, 2, 3, 5, 8))
typeInfo(myDF)
myDF <- myDF[-1, ]
myDF

myDF2 <- as.numeric(myDF)
myDF2     # Whoa! what just happened ?

myDF3 <- as.numeric(as.character(myDF))
myDF3     # :-)
~~~

 Task:
 Repair the sup3 data.frame - realod it with stringsAsFactors = FALSE

~~~r
sup3[1:10,]
head(sup3)
tail(sup3)
nrow(sup3)
ncol(sup3)
sup3$genes[1:10]
sup3$genes[1:10]
~~~



 Now we can finally return to our original question and try e.g ...

~~~r
sup3[sup3$genes == "Cd19", ]
sup3[sup3$genes == "Lyz2", ]

# But!
sup3[sup3$genes == "B220", ]    Not found!
~~~

 B220 has a number of synonyms as a quick Google search shows:
B220
CD45
CD45R
GP180
L-CA
LCA
LY5
PTPRC
T200

 Is there a convenient way to convert such lines into a character vector?
 Yes: as we have done above, define the list as one string constant, then use strsplit() on it.

~~~r
s <- "B220
CD45
CD45R
GP180
L-CA
LCA
LY5
PTPRC
T200"
B220synonyms <- unlist(strsplit(s, "\n"))

?"%in%"
B220synonyms %in% toupper(sup3$genes)
~~~

 However - no luck. None of the synonyms is in the table either.
 How do we know that this is not a problem with our expression?

 Positive control!
 
~~~r
c(B220synonyms, "CD19") %in% toupper(sup3$genes)
~~~

 Task: check if our "characteristic genes" are all in the table then find the enrichment vectors for the subset Bst2, Siglech, Ly6d, Irf8


 To do more than that, we really need to look at writing "programs".

