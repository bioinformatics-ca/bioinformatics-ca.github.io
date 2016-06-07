# ==============================================================================
# Exploratory Analysis of Biological Data using R, 2015
# Integrated Assignment Part 1 Answers
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
# SECTION 1  Set your working directory
#
# e.g. (BUT MAKE IT YOUR OWN!)
setwd("C:/Users/Owner/Desktop/Goldenberg Lab/OICR Workshop Materials/integratedassignmentfromlastyear/")

# ==============================================================================
# SECTION 2  Import the data into R
#
load("ccleCgc.rda")
ls()

# ==============================================================================
# SECTION 3  Examine the environment and objects
# 

# 3a  Get the list of objects in the workspace using ls().
ls()

# 3b  What is the difference between ls() and list.files()?
list.files()
# Answer:  ls() lists objects in the R environment and list.files() lists the
#          files stored on the disk.

# For now, we will only work with `expr`, `cn`, and `pheno`.

# 3c  What are the classes of `expr`, `cn`, and `pheno`?
#     e.g. class(expr) outputs the class of `expr`
class(expr)
class(cn)
class(pheno)

# 3d  How many columns and rows does `expr` have? (Hint: nrow(), ncol())
nrow(expr)
ncol(expr)
# Alternatively
dim(expr)

# 3e  How many columns and rows do the data frames `cn`, `expr`, `pheno` have?
dim(cn)
dim(expr)
dim(pheno)

# 3f  What are the names of the rows of `expr`? The names of the columns?
rownames(expr)
colnames(expr)

# 3g  What are the names of the rows of `pheno`? The names of the columns?
rownames(pheno)
colnames(pheno)


# ==============================================================================
# SECTION 4  Convert and rearrange the data
# 

# 4a  Restructure the data into an easy-to-use format.
# Already done.

# 4b  Examine the description of the data by typing `?ccleCgc`.
#     Type `q` to get back to the command line, after reading the description.
#     For all the data frames (e.g. `cn`), what do the rows represent?
# Answer:  Each row represents a cancer cell line.

# We should make sure that the data frames describe the same cell lines
# and that the data are arranged in the same order.
# We can do so by ensuring that the row names are the same.
# The code below does this check for for the data frames `cn` and `expr`.
identical(rownames(cn), rownames(expr))

# 4c  Now repeat this comparison between `cn` and `pheno`.
identical(rownames(cn), rownames(pheno))


# ==============================================================================
# SECTION 5  Examine the data distributions
# 

# 5a  Examine the usage instructions for hist(), using ? hist.
# Note:  hist(x) does not work if x is a data.frame.
#        x must be a vector. Matrices are automatically converted to vectors.
#        as.matrix(x) converts x into a matrix.
#     Plot a histogram of `expr`.
hist(as.matrix(expr))

# 5b  Plot a prettier histogram of `expr`.
#     Use 50 breaks and set the colour of the bars to "skyblue".
#     Additionally, set the pararmeter freq to FALSE.
# Hint:  In the usage instructions of hist(), what do the parameters col, 
#        breaks, and freq do?
# Note:  Don't forget to convert `expr` into a matrix.
hist(as.matrix(expr), breaks=50, col="dodgerblue", freq=FALSE)

# If hist() have been called with freq=FALSE, you can superimpose a 
# density curve by the following command:
lines(density(as.matrix(expr)), col="red",lwd=3)

# 5c  Plot a histogram of `cn` with 100 breaks and orange bars.
#     Superimpose a density plot.
hist(as.matrix(cn), breaks=100, freq=FALSE, col="darkorchid")
lines(density(as.matrix(cn)), col="red",lwd=3)

# 5d  Create a quantile-quantile plot of `cn` against the normal 
#     distribution. Draw a reference line through the data.
# Hint:  Use qqnorm() to compare against a normal distribution.
qqnorm(as.matrix(cn), pch='.')
qqline(as.matrix(cn), col="red")

# You can also superimpose a normal distribution on top of a histogram:
z <- as.numeric(as.matrix(cn))
hist(z, breaks=100, freq=FALSE, col="deeppink2")
curve(dnorm(x, mean=mean(z), sd=sd(z)), col="blue", add=TRUE,lwd=3)


# ==============================================================================
# SECTION 6  Explore the gene expression data
# 

# 6a  Optional: Check if the tumor suppressor gene TP53 is in the dataset. 
# Hint:  There are two ways: one using `%in%` and the other using `which()`.
#        Recall whether the genes are along rows or columns of the data.
"TP53" %in% colnames(expr)
which(colnames(expr)=="TP53")

# 6b  Print the expression values of TP53 across all cell lines.
expr[,"TP53"]

# 6c  What are the main statistics?
#     (e.g. min, max, median, mean, standard deviation)
summary(expr[,"TP53"])
# looking at them one-by-one:
min(expr[,"TP53"])
max(expr[,"TP53"])
median(expr[,"TP53"])
mean(expr[,"TP53"])
sd(expr[,"TP53"])

# 6d  Look at the histogram and QQ-plot of TP53 expression values.
#     Are the values normally distributed?
hist(expr[,"TP53"])
qqnorm(expr[,"TP53"])
qqline(expr[,"TP53"], col="red",lwd=3)

# 6e  From the histogram, we can see that some cell lines have low TP53 
#     expression, while others have high TP53 expression.
#     Let's consider expression values of 5.5 or lower as low expression.
#     First, check which entries of the TP53 expression values are < 5.5.
expr[,"TP53"] < 5.5
#    How many cell lines with low TP53 expression are there? 
#    Use the answer above in place of the ???? in the below partial code:
nrow(expr[expr[,"TP53"] < 5.5,])    # nrow(expr[????,])

# 6f  What are the names of these cell lines?
rownames(expr[expr[,"TP53"] < 5.5,])

# 6g  Create a scatterplot of TP53 copy-number and expression.
#     What do you notice about the trend in the values?  TP53 encodes a protein
#     involved in DNA damage response TP53 is a tumour suppressor gene because
#     it induces programmed cell death or cell cycle arrest upon DNA damage.
  ## try manipulating the point type and color
  
plot(cn[,"TP53"], expr[,"TP53"],
     pch=19,col='forestgreen')

# Let's make a boxplot of TP53 expression value across different cancer sites.

# First, make a data frame of the site and T53 expression value 
tp53.df <- data.frame(
  site = pheno$site,
  expr = expr[,"TP53"]
)

# 6h  Print out the first 10 rows of this data frame.
#     Notice that some of the cell lines come from the same site (e.g. skin).
tp53.df[1:10,]

# Then, adjust the margins of the plot and make the labels horizontal
par(mar = c(16, 4, 1, 2), las=2)

# Finally, create box plot, grouping by site
boxplot(tp53.df$expr ~ tp53.df$site)

# ==============================================================================
# SECTION 6.5  Creating publication quality plots
# 

# Import the ggplot2 library
library(ggplot2)
library(ggthemes)
# 6.5a  Create a data frame containing TP53 and MYC expression data
#      with each row representing one sample.
d <- expr[, c("TP53", "MYC")]

# 6.5b  Create an ggplot object and specify plot variables.
(g <- ggplot(d, aes(x=TP53, y=MYC)) + geom_point())
print(g)

# 6.5c  Add plot title, change point size and color 
#       and theme specification to the ggplot object.
(g <- ggplot(d, aes(x=TP53, y=MYC)) + 
   geom_point(col="dodgerblue3",size=4) + 
   theme_economist_white() + ggtitle("TP53 vs MYC Cell line Expression Levels"))

# ==============================================================================
# SECTION 7  Principal component analysis
# 

# 7a  Apply PCA on the expression data.
# Note:  prcomp() expects the samples to be along the rows
expr.pr <- prcomp(expr)

# 7b  Plot the first two principal components.
plot(expr.pr$x[,1], expr.pr$x[,2],
     col='firebrick',pch=19)

# 7c  What are the groups? Do they represent different cancer types?
# Hint: You can plot data points in different colours by specifiying the `col`
#       parameter using a vector of numbers.
#       The phenotype information for the cell lines are in `pheno`.
plot(expr.pr$x[,1], expr.pr$x[,2], col=as.numeric(pheno$site),
	xlim=c(-10, 90))

# 7d  Add a legend to the plot.
sites <- levels(pheno$site)
legend("topright", legend=sites, fill=1:length(sites), bty="n")

# 7e  There are too many different cancer types.
#     Let's subset a few: breast, prostate, ovary, lung, skin, bone, 
#     haematopoietic_and_lymphoid_tissue, central_nervous_system.
#     Create a vector of cell lines correpsonding to these cancer types.
cell.lines <- rownames(pheno)[pheno$site %in%
	c("breast", "prostate", "ovary", "lung", "skin", "bone",
	"haematopoietic_and_lymphoid_tissue", "central_nervous_system")]

# 7f  Create a subset of `pheno` and call it `pheno.sub`.
pheno.sub <- pheno[cell.lines, ]

# Note:  You need to re-create factor variables to remove missing factor levels.
pheno.sub$site <- factor(pheno.sub$site)

# 7g  Create a subset of `expr` and apply PCA.
expr.sub <- expr[cell.lines, ]
expr.sub.pr <- prcomp(expr.sub)

# 7h  Plot the first two principal components of the data subset.
plot(expr.sub.pr$x[,1], expr.sub.pr$x[,2], col=as.numeric(pheno.sub$site),
     xlim=c(-30, 15), ylim=c(-20, 15))
sites <- levels(pheno.sub$site)
legend("bottomleft", legend=sites, fill=1:length(sites), bty="n")



par(mfrow=c(1,1))
plot(expr.sub.pr$x[,1], expr.sub.pr$x[,2], 
     col=c("blue","red","forestgreen","darkorchid3","deeppink1","dodgerblue","black","orange")[factor(as.numeric(pheno.sub$site))],
	xlim=c(-30, 15), ylim=c(-20, 15),pch=19)
sites <- levels(pheno.sub$site)
legend("bottomleft", legend=sites, fill=c("blue","red","forestgreen","darkorchid3","deeppink1","dodgerblue","black","orange"), bty="n")


# ==============================================================================
# SECTION 8  Explore the gene expression again
# 

# 8a  Compute the correlation coefficient between the expression values of
#     RUNX1 and JUN using cor().
cor(expr[,"RUNX1"], expr[,"JUN"])

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
  values1 <- expressionMatrix[,gene1]
  values2 <- expressionMatrix[,gene2]
  cor(values1, values2)
}

# You will now build a FOR-loop as shown in the lecture notes to correlate 
# RUNX1 expression values to expression of other genes.

# 8c  Which genes do you want to correlate with RUNX1?
#          Assign that to a vector called testedGenes. 
testedGenes <- colnames(expr)

# 8d  How many genes do we have?
#     Store this value in the variable called numberOfGenes.
numberOfGenes <- length(testedGenes)

# 8e  Create an vector of zeroes to store correlation results (corResults).
#     The vector should contain one element for every gene tested. 
corResults <- rep(0, length(testedGenes))

# 8f  Create a loop that calls the function correlateTwoGenes
#     consecutively on all genes in the matrix. 
#     The first argument gene1 of correlateTwoGenes is always "RUNX1",
#     but the second argument gene2 is given by the for-loop from the
#     vector testedGenes. The third argument of the function is always the
#     expression matrix (expr). 
#     Assign the computed correlation coefficient to the corresponding 
#     value in the corResults vector. 
#     The for-loop should iterate from 1 to numberOfGenes. 
for (i in 1:numberOfGenes) {
  corResults[i] <- correlateTwoGenes("RUNX1", testedGenes[i], expr)
}

# 8g  Select the genes that have a moderate positive correlation with 
#     RUNX1 (correlation of 0.3 more). 
testedGenes[corResults>=0.3]

# You should see six genes (TAL1, CCND3, CHCHD7, ETV6, LYL1, and RUNX1).
# One of the top correlated genes is ETV6, another transcription factor
# involved in development and cancer. RUNX1 and ETV6 are known to make up
# fusion proteins through chromosomal translocations, leading to leukemia.

# 8h  Investigate the correlation of RUNX1 and ETV6.
#     Plot their expression values as a scatterplot 
#     (i.e. expression of ETV6 vs. expression of RUNX1)
plot(expr[,"RUNX1"], expr[,"ETV6"])

# 8i  Compute correlation p-value using cor.test().
#     What is the p-value of this association?
h <- cor.test(expr[,"RUNX1"], expr[,"ETV6"])
h$p.value
h$estimate

# 8j  Create a linear regression model to predict ETV6 values
#     from RUNX1 values. 
# Note:  The arguments of plot(x, y) are reversed in lm(y ~ x). 
#        So RUNX1 and ETV6 should switch places.
etvModel <- lm(expr[,"ETV6"] ~ expr[,"RUNX1"])

# 8k  Investigate the regression model using the summary() function.
#     Compare the gained p-value with the p-value from cor.test.
summary(etvModel)

# 8l  Add the linear regression line on the plot using a thick red line.
abline(etvModel, col="red", lwd=3)


# ==============================================================================
# SECTION 9  Explore the pharmacological profiles
# 

# 9a  Four response variables were measured (ic50, ec50, act.area, and act.max).
#     These data are contained with the `pharm` list.
#     Create scatter plots of each pair of response variables.
plot(log10(as.matrix(pharm$ic50)), as.matrix(pharm$act.max), pch=20)
plot(log10(as.matrix(pharm$ec50)), as.matrix(pharm$act.max), pch=20)
plot(log10(as.matrix(pharm$ec50)), as.matrix(pharm$act.area), pch=20)

# 9b  Can you see different clusters in the act.max vs ec50 plot?
#     What do each cluster represent?
# Hint:  Zoom into the region where most data points lie using xlim and ylim
plot(log10(as.matrix(pharm$ec50)), as.matrix(pharm$act.max),
     pch=20, xlim=c(-3, 1), ylim=c(-100, 0))

# 9c  Paclitaxel is a non-targeted anti-cancer agent, whereas PLX4720 is a 
#     targeted anti-cancer agent (kinase inhibitor designed to inhibit the
#     oncogene BRAF^V600E).
#     How do Paclitaxel and PLx4720 compare in terms of their activities 
#     across cell lines?
#     Plot a scatter plot for Paclitaxel. Do the same for PLX4720.
# Note:  We've been provided with a 'plot.R' file with a 
# plot.drug.activity function within it, first 'source()' the plot.R function
source('plot.R')
plot.drug.activity(pharm, "Paclitaxel")
plot.drug.activity(pharm, "PLX4720")

# 9d  Repeat for other drugs.
#     Panobinostat (HDAC inhibitor)
plot.drug.activity(pharm, "Panobinostat")
#     17-AAG (Antibiotic drug re-purposed for fighting cancer)
plot.drug.activity(pharm, "X17.AAG")
#     Nutlin-3 (Drug designed to restore wildtype p53 function)
plot.drug.activity(pharm, "Nutlin.3")
#     Erlotinib (Kinase inhibitor designed to inhibit EGFR)
plot.drug.activity(pharm, "Erlotinib")

# 9e  The PLX4720 drug is designed to inhibit the cancers harbouring 
#     mutant BRAF with the V600E substitution. Base on the pharmacological
#     profiles of different cell lines, does this drug appear effective
#     against its intended target?
# Hint:  Compare the response of cell lines with BRAF^V600E to those without.
#        Assign colours to each point based on BRAF^V600E mutation status.
plot.drug.activity(pharm, "PLX4720", group=mut$BRAF.V600E)
# Answer:  Cell lines with BRAF^V600E mutation are more susceptible to PLX4720.

# 9f  Are the chemotherapeutic agents, Paclitaxel, Irinotecan, and Topotecan,
#     universally effective against all cancer types?
plot.drug.activity(pharm, "Paclitaxel", group=pheno$site)
plot.drug.activity(pharm, "Irinotecan", group=pheno$site)
plot.drug.activity(pharm, "Topotecan", group=pheno$site)


# ==============================================================================
# SECTION 10  Create publication-quality plots 
# 

# Recall the plot 'g' we created above. We can save the plot to a file using qdraw() from the 
# io package, which will save files in the correct formats as specified by
# the file extensions in the file names. Additionally, qdraw() creates smaller
# drawing areas by default to ensure the labels are relatively larger.
# The dimensions of the plot can be easily set by changing the width and
# height parameters of qdraw().
library('io')
qdraw(g, file="myc-tp53.pdf")
qdraw(g, file="myc-tp53.png")

# 10a  Remake the plot and show the histotypes in different colours.
d <- cbind(expr[, c("TP53", "MYC")], histotype=pheno$histotype)
g <- ggplot(d, aes(x=TP53, y=MYC, colour=histotype)) + 
  geom_point(size=4) + theme_bw()
qdraw(g, width=9)

# 10b  Show the same plot but only for cancers of haematopoietic or
#      lymphoid origin
d <- cbind(expr[, c("TP53", "MYC")], pheno[, c("site", "histotype")])
g <- ggplot(d[d$site == "haematopoietic_and_lymphoid_tissue", ],
            aes(x=TP53, y=MYC, colour=histotype)) + geom_point(size=4) + theme_bw()
qdraw(g, width=7)

# 10c  Visualize the expression data by PCA.
#      Apply PCA to the entire dataset.
#      Show only the cancer types breast, skin, lung, and liver.
expr.pr <- prcomp(expr)
d <- cbind(expr.pr$x, pheno[, c("site", "histotype")])
g <- ggplot(d[d$site %in% c("breast", "skin", "lung", "liver"), ],
            aes(x=PC1, y=PC2, colour=site)) + geom_point(size=4) + theme_bw()
qdraw(g)

# 10d  Repeat the same PCA analysis with copy-number data.
cn.pr <- prcomp(cn)
d <- cbind(cn.pr$x, pheno[, c("site", "histotype")])
g <- ggplot(d[d$site %in% c("breast", "skin", "lung", "liver"), ], 
            aes(x=PC1, y=PC2, colour=site)) + geom_point(size=4) + theme_bw()
qdraw(g)
