# ==============================================================================
# Exploratory Analysis of Biological Data using R, 2015
# Integrated Assignment Part 1 Questions
#
# @Authors:  David JH Shih <djh.shih@gmail.com>
#            Catalina Anghel <catalina.anghel@oicr.on.ca>
# @License:  GNU General Public License v3 
# @Created:  2015-05-11
#
# @Input:    Data from the Cancer Cell Line Encyclopedia (CCLE)
#            Downloaded from: http://www.broadinstitute.org/ccle
#            Associated publication:
#
#            Barretina et al. The Cancer Cell Line Encylopedia enables
#            predictive modelling of anticancer drug sensitivity.
#            Nature 483, 603-607 (2012).
#
# @Output:   Statistics, Plots


# ==============================================================================
# SECTION 1  Install the required R package
#

# You should have built and installed the CLLE package.


# ==============================================================================
# SECTION 2  Import the data into R
#

# 2a  Load the "CCLE" R library by using the library() function.
setwd("YOUR WORKING DIRECTORY")

# 2b  Import the data into the workspace
load("ccleCgc.rda");


# ==============================================================================
# SECTION 3  Examine the environment and objects
# 

# 3a  Get the list of objects in the workspace using ls().
????

# 3b  What is the difference between ls() and list.files()?
????

# For now, we will only work with `expr`, `cn`, and `pheno`.

# 3c  What are the classes of `expr`, `cn`, and `pheno`?
#     e.g. class(expr) outputs the class of `expr`
????
????
????

# 3d  How many columns and rows does `expr` have? (Hint: nrow(), ncol())
????
????

# 3e  How many columns and rows do the data frames `cn`, `expr`, `pheno` have?
????
????
????

# 3f  What are the names of the rows of `expr`? The names of the columns?
????
????

# 3g  What are the names of the rows of `pheno`? The names of the columns?
????
????


# ==============================================================================
# SECTION 4  Convert and rearrange the data
# 

# 4a  Restructure the data into an easy-to-use format.
# *** Already done. ****

# 4b  Each row of the data represents a cancer cell line.

# We should make sure that the data frames describe the same cell lines
# and that the data are arranged in the same order.
# We can do so by ensuring that the row names are the same.
# The code below does this check for for the data frames `cn` and `expr`.
identical(rownames(cn), rownames(expr));

# 4c  Now repeat this comparison between `cn` and `pheno`.
????

# ==============================================================================
# SECTION 5  Examine the data distributions
# 

# 5a  Examine the usage instructions for hist(), using ? hist.
# Note:  hist(x) does not work if x is a data.frame.
#        x must be a vector. Matrices are automatically converted to vectors.
#        as.matrix(x) converts x into a matrix.
#     Plot a histogram of `expr`.
????

# 5b  Plot a prettier histogram of `expr`.
#     Use 50 breaks and set the colour of the bars to "skyblue".
#     Additionally, set the pararmeter freq to FALSE.
# Hint:  In the usage instructions of hist(), what do the parameters col, 
#        breaks, and freq do?
# Note:  Don't forget to convert `expr` into a matrix.
????

# If hist() have been called with freq=FALSE, you can superimpose a 
# density curve by the following command:
lines(density(as.matrix(expr)), col="red");

# 5c  Plot a histogram of `cn` with 100 breaks and orange bars.
#     Superimpose a density plot.
????
????

# 5d  Create a quantile-quantile plot of `cn` against the normal 
#     distribution. Draw a reference line through the data.
# Hint:  Use qqnorm() to compare against a normal distribution.
????
????

# You can also superimpose a normal distribution on top of a histogram:
z <- as.numeric(as.matrix(cn));
hist(z, breaks=100, freq=FALSE, col="orange");
curve(dnorm(x, mean=mean(z), sd=sd(z)), col="blue", add=TRUE);


# ==============================================================================
# SECTION 6  Explore the gene expression data
# 

# 6a  Optional: Check if the tumor suppressor gene TP53 is in the dataset. 
# Hint:  There are two ways: one using `%in%` and the other using `which()`.
#        Recall whether the genes are along rows or columns of the data.
????

# 6b  Print the expression values of TP53 across all cell lines.
????

# 6c  What are the main statistics?
#     (e.g. min, max, median, mean, standard deviation)
????

# 6d  Look at the histogram and QQ-plot of TP53 expression values.
#     Are the values normally distributed?
????

# 6e  From the histogram, we can see that some cell lines have low TP53 
#     expression, while others have high TP53 expression.
#     Let's consider expression values of 5.5 or lower as low expression.
#     First, check which entries of the TP53 expression values are < 5.5.
????
#    How many cell lines with low TP53 expression are there? 
#    Use the answer above in place of the ???? in the below partial code:
nrow(expr[????,])

# 6f  What are the names of these cell lines?
????

# 6g  Create a scatterplot of TP53 copy-number and expression.
#     What do you notice about the trend in the values?  TP53 encodes a protein
#     involved in DNA damage response; TP53 is a tumour suppressor gene because
#     it induces programmed cell death or cell cycle arrest upon DNA damage.
????

# Let's make a boxplot of TP53 expression value across different cancer sites.

# First, make a data frame of the site and T53 expression value 
tp53.df <- data.frame(
  site = pheno$site,
  expr = expr[,"TP53"]
);

# 6h  Print out the first 10 rows of this data frame.
#     Notice that some of the cell lines come from the same site (e.g. skin).
????

# Then, adjust the margins of the plot and make the labels horizontal
par(mar = c(16, 4, 1, 2), las=2);

# Finally, create box plot, grouping by site
boxplot(tp53.df$expr ~ tp53.df$site);


# ==============================================================================
# SECTION 6.5  Creating publication quality plots
# 
## First library in ggplot2
library('ggplot2')

# 6.5a  Create a data frame containing TP53 and MYC expression data
#      with each row representing one sample.
d <- expr[, c("TP53", "MYC")];

# 6.5b  Create an ggplot object and specify plot variables.
g <- ggplot(d, aes(x=TP53, y=MYC)) + geom_point();
????
????

# 6.5c  Add plot title, change point size and color
#       and theme specification to the ggplot object.
library(ggthemes)
???
???


# ==============================================================================
# SECTION 7  Principal component analysis
# 

# 7a  Apply PCA on the expression data.
# Note:  prcomp() expects the samples to be along the rows
????

# 7b  Plot the first two principal components.
plot(expr.pr$x[,????], expr.pr$x[,????],

# 7c  What are the groups? Do they represent different cancer types?
# Hint: You can plot data points in different colours by specifiying the `col`
#       parameter using a vector of numbers.
#       The phenotype information for the cell lines are in `pheno`.
plot(????, ????, col=as.numeric(pheno$site),
	xlim=c(-10, 90));

# 7d  Add a legend to the plot.
sites <- levels(pheno$site);
legend("topright", legend=????, fill=1:length(????), bty="n");

# 7e  There are too many different cancer types.
#     Let's subset a few: breast, prostate, ovary, lung, skin, bone, 
#     haematopoietic_and_lymphoid_tissue, central_nervous_system.
#     Create a vector of cell lines correpsonding to these cancer types.
cell.lines <- rownames(pheno)[????];

# 7f  Create a subset of `pheno` and call it `pheno.sub`.
pheno.sub <- ????;

# Note:  You need to re-create factor variables to remove missing factor levels.
pheno.sub$site <- factor(pheno.sub$site);

# 7g  Create a subset of `expr` and apply PCA.
expr.sub <- ????;
expr.sub.pr <- ????;

# 7h  Plot the first two principal components of the data subset.
plot(????, ????, col=????,
	xlim=c(-30, 15), ylim=c(-20, 15));
sites <- levels(pheno.sub$site);
legend("bottomleft", legend=????, fill=????, bty="n");


# ==============================================================================
# SECTION 8  Explore the gene expression again
# 

# 8a  Compute the correlation coefficient between the expression values of
#     RUNX1 and JUN using cor().
????

# The expressions of RUNX1 and JUN are not very correlated.
# Let's automate the procedure to test correlations between RUNX1 and
# other genes to find better candidates.

# 8b  Write a function correlateTwoGenes(gene1, gene2, expressionMatrix)
#     to find the correlation between two given genes in our dataset. 
#     The function has three parameters: gene1, gene2, expressionMatrix.
#     The code in the function needs to extract gene-specific expressions
#     from the matrix and compute their correlation.
#     The output of the function is the Pearson correlation coefficient of
#     the expression values of the two genes. 
# Example:  correlateTwoGenes("RUNX1", "TP53", expr) should give: 0.07595393
correlateTwoGenes <- function(gene1, gene2, expressionMatrix) {
  values1 <- ????;
  values2 <- ????;
  ????
}

# You will now build a FOR-loop as shown in the lecture notes to correlate 
# RUNX1 expression values to expression of other genes.

# 8c  Which genes do you want to correlate with RUNX1?
#          Assign that to a vector called testedGenes. 
????

# 8d  How many genes do we have?
#     Store this value in the variable called numberOfGenes.
????

# 8e  Create an vector of zeroes to store correlation results (corResults).
#     The vector should contain one element for every gene tested. 
corResults <- rep(0, length(????));

# 8f  Create a loop that calls the function correlateTwoGenes
#     consecutively on all genes in the matrix. 
#     The first argument gene1 of correlateTwoGenes is always "RUNX1",
#     but the second argument gene2 is given by the for-loop from the
#     vector testedGenes. The third argument of the function is always the
#     expression matrix (expr). 
#     Assign the computed correlation coefficient to the corresponding 
#     value in the corResults vector. 
#     The for-loop should iterate from 1 to numberOfGenes. 
for (????) {
  ????
}

# 8g  Select the genes that have a moderate positive correlation with 
#     RUNX1 (correlation of 0.3 more). 
????

# You should see six genes (TAL1, CCND3, CHCHD7, ETV6, LYL1, and RUNX1).
# One of the top correlated genes is ETV6, another transcription factor
# involved in development and cancer. RUNX1 and ETV6 are known to make up
# fusion proteins through chromosomal translocations, leading to leukemia.

# 8h  Investigate the correlation of RUNX1 and ETV6.
#     Plot their expression values as a scatterplot 
#     (i.e. expression of ETV6 vs. expression of RUNX1)
????

# 8i  Compute correlation p-value using cor.test().
#     What is the p-value of this association?
????

# 8j  Create a linear regression model to predict ETV6 values
#     from RUNX1 values. 
# Note:  The arguments of plot(x, y) are reversed in lm(y ~ x). 
#        So RUNX1 and ETV6 should switch places.
????

# 8k  Investigate the regression model using the summary() function.
#     Compare the gained p-value with the p-value from cor.test.
????

# 8l  Add the linear regression line on the plot using a thick red line.
????


# ==============================================================================
# SECTION 9  Explore the pharmacological profiles
# 

# 9a  Four response variables were measured (ic50, ec50, act.area, and act.max).
#     These data are contained with the `pharm` list.
#     Create scatter plots of each pair of response variables.
????
????
????

# 9b  Can you see different clusters in the act.max vs ec50 plot?
#     What do each cluster represent?
# Hint:  Zoom into the region where most data points lie using xlim and ylim
????

# 9c  Paclitaxel is a non-targeted anti-cancer agent, whereas PLX4720 is a 
#     targeted anti-cancer agent (kinase inhibitor designed to inhibit the
#     oncogene BRAF^V600E).
#     How do Paclitaxel and PLx4720 compare in terms of their activities 
#     across cell lines?
#     Plot a scatter plot for Paclitaxel. Do the same for PLX4720.
# Note:  Within the directory we've saved a script with a specialized plot function
# called "plot.drug.activity". Use the command source() to read this function
# into your environment and use it to plot the drug activities.
source("plot.R")

????
????

# 9d  Repeat for other drugs.
#     Panobinostat (HDAC inhibitor)
????
#     17-AAG (Antibiotic drug re-purposed for fighting cancer)
????
#     Nutlin-3 (Drug designed to restore wildtype p53 function)
????
#     Erlotinib (Kinase inhibitor designed to inhibit EGFR)
????

# 9e  The PLX4720 drug is designed to inhibit the cancers harbouring 
#     mutant BRAF with the V600E substitution. Base on the pharmacological
#     profiles of different cell lines, does this drug appear effective
#     against its intended target?
# Hint:  Compare the response of cell lines with BRAF^V600E to those without.
#        Assign colours to each point based on BRAF^V600E mutation status.
????
# Answer:  Cell lines with BRAF^V600E mutation are more susceptible to PLX4720.

# 9f  Are the chemotherapeutic agents, Paclitaxel, Irinotecan, and Topotecan,
#     universally effective against all cancer types?
????
????
????

# ==============================================================================
# SECTION 10  Create MORE publication-quality plots 
# 

# Import the ggplot2 library
library(ggplot2);


# Alternatively, we can save the plot to a file using qdraw() from the 
# io package, which will save files in the correct formats as specified by
# the file extensions in the file names. Additionally, qdraw() creates smaller
# drawing areas by default to ensure the labels are relatively larger.
# The dimensions of the plot can be easily set by changing the width and
# height parameters of qdraw().
library(io);
qdraw(g, file="myc-tp53.pdf");
qdraw(g, file="myc-tp53.png");

# 10e  Remake the plot and show the histotypes in different colours.
d <- cbind(expr[, c("TP53", "MYC")], histotype=pheno$histotype);
????
????
????

# 10f  Show the same plot but only for cancers of haematopoietic or
#      lymphoid origin
d <- cbind(expr[, c("TP53", "MYC")], pheno[, c("site", "histotype")]);
????
????
????

# 10g  Visualize the expression data by PCA.
#      Apply PCA to the entire dataset.
#      Show only the cancer types breast, skin, lung, and liver.
????
????
????
????
????

# 10h  Repeat the same PCA analysis with copy-number data.
????
????
????
????
????
