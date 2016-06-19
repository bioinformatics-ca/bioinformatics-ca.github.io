## R review code
# Modified by Florence Cavalli (2013, 2014, 2016), and initially written by Richard De Borja and Cindy Yao (2012) with some additions by Sorana Morrissy (2015)

### SLIDE 7 ######################################################################################
2+5
30/10
20-5
log10(20)
pi
exp(-2)
sin(0.4)

### SLIDE 8 ######################################################################################
weight.a <- 10
weight.a
weight.b <- 30
weight.b
total.weight <- weight.a + weight.b
total.weight

### SLIDE 10 ######################################################################################
# get the current working directory:
getwd()
# set a new working directory:
setwd("C:/myPATH")
setwd("~/myPATH") # on Mac
setwd("/Users/florence/myPATH") # on Mac
# list the files in the current working directory:
list.files()
# list objects in the current R session:
ls()


### SLIDE 11 ######################################################################################
help(sum)
# q to quit the help page
?sum

help.search("plot")
??plot

### SLIDE 16 ######################################################################################
## vectors
numeric.vector <- c(1,2,3,4,5,6,2,1)
numeric.vector

character.vector <- c("Fred", "Barney", "Wilma", "Betty")
character.vector

logical.vector <- c(TRUE, TRUE, FALSE, TRUE)
logical.vector


### SLIDE 17 ######################################################################################
character.vector
character.vector[2]
character.vector[2:3]
character.vector[c(2,4)]

### SLIDE 18 ######################################################################################
# check the class of the object
class(numeric.vector)
class(character.vector)
class(logical.vector)

# check the structure of the object
str(logical.vector)

### SLIDE 19 ######################################################################################
gender <- c(1,2,1,1,1,2)
gender

gender.factor <- as.factor(gender)
gender.factor

levels(gender.factor) <- c("male", "female")
gender.factor

### SLIDE 20 ######################################################################################
gender.factor
gender.factor[2]
gender.factor[2:4]
gender.factor[c(1,4)]

### SLIDE 21 ######################################################################################
## matrix
?matrix
matrix.example <- matrix(1:12, nrow = 3, ncol=4, byrow = FALSE)
matrix.example

matrix.example <- matrix(1:12, nrow = 3, ncol=4, byrow = TRUE)
matrix.example


### SLIDE 22-23 ######################################################################################
dataset.a <- c(1,22,3,4,5)
dataset.b <- c(10,11,13,14,15)
dataset.a
dataset.b

rbind.together <- rbind(dataset.a, dataset.b)
rbind.together

cbind.together <- cbind(dataset.a, dataset.b)
cbind.together

### SLIDE 24 ######################################################################################
matrix.example
matrix.example[2,4]
matrix.example[2,]
matrix.example[,4]

### SLIDE 25 ######################################################################################
colnames(matrix.example) <- c("Sample1","Sample2","Sample3","Sample4")
rownames(matrix.example) <- paste("gene",1:3,sep="_")
matrix.example

matrix.example[,"Sample2"]
matrix.example[1,"Sample2"]
matrix.example["gene_1","Sample2"]

### SLIDE 26 ######################################################################################
## data.frame
people.summary <- data.frame(
                             age = c(30,29,25,25),
                             names = c("Fred", "Barney", "Wilma", "Betty"),
                             gender = c("m", "m", "f", "f")
                             )
people.summary

### SLIDE 27 ######################################################################################
people.summary
people.summary[2,1]
people.summary[2,]
people.summary[,1]
people.summary$age

### SLIDE 28-29 ######################################################################################
## List
together.list <- list(
                      vector.example = dataset.a, 
                      matrix.example = matrix.example,
                      data.frame.example = people.summary
                      )
together.list

### SLIDE 30 ######################################################################################
together.list$matrix.example
together.list$matrix.example[,3]
together.list["matrix.example"]
together.list[["matrix.example"]]
together.list[["matrix.example"]][,2]


### SLIDE 31 ######################################################################################
## Basic functions
sum(c(1,2,3))
log2(10)
sin(0.24)
mean(c(1,2,3,4,5))

### SLIDE 32 ######################################################################################
character.vector
length(character.vector)
matrix.example
nrow(matrix.example)
ncol(matrix.example)
dim(matrix.example)

### SLIDE 33 ######################################################################################
?read.table
## read.table("myDataFile.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE) ##Example how to read a text file, not ran here
?write.table
write.table(people.summary, file="File_name_people_summary.txt", quote=FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)


### solution Exercises I - SLIDE 34 ###############################################################
### SLIDE 35 ######################################################################################
# Use the assignment operator to store three values of your choice.
value.a <- 10
value.b <- 3
value.c <- 12

# Calculate and store the sum of the values from above.
sum.values <- value.a + value.b + value.c
sum.values

### Solution Exercises II - SLIDE 36 #############################################################
### SLIDE 37 ######################################################################################
# What data type is the cars dataset? and its dimensions?
class(cars)
dim(cars)
## or
str(cars)

## To look at the object
head(cars) ## print a subset of the dataset

## Access to the speed values
cars$speed

### SLIDE 38 ######################################################################################
# Access and store only the cars data with speeds greater than 15. How many cars does this affect?
cars.greater.speed <- cars[cars$speed > 15,]
nrow(cars.greater.speed)

# Reformat the cars data into a list
cars.as.list <- list(
	SPEED = cars$speed,
	DISTANCE = cars$dist
	)
cars.as.list
names(cars.as.list)

### SLIDE 39 ######################################################################################
# Access only the cars data with speeds greater than 15 from the list you just created. How many cars does this affect? Did you get the same results as above?
cars.as.list.greater.speed <- cars.as.list$SPEED[cars.as.list$SPEED > 15]
cars.as.list.greater.speed
length(cars.as.list.greater.speed)

### SLIDE 40 ######################################################################################
# What does stringsAsFactors in the data.frame function do?
?data.frame

### SLIDE 42 ######################################################################################
plot(x=cars$speed, y=cars$dist)

plot(x=cars$speed, y=cars$dist,
    xlab = "Speed",
    ylab = "Distance",
    cex.lab = 1.5,
    main = "A nice scatter plot",
    pch = 16,
    bty = "n",
    col = "dark blue",
    las = 1
    )
## las
## How to change the axes label style in R
## To change the axes label style, use the graphics option las (label style). This changes the orientation angle of the labels:
## 0: The default, parallel to the axis
## 1: Always horizontal
## 2: Perpendicular to the axis
## 3: Always vertical

## bty
## To change the type of box round the plot area, use the option bty (box type):
## "o": The default value draws a complete rectangle around the plot.
## "n": Draws nothing around the plot.

### SLIDE 43 ######################################################################################
hist(cars$speed)

hist(cars$speed,
    xlab = "Speed",
    ylab = "Number of cars",
    cex.lab = 1.5,
    main = "A nice histogram",
    col = "cyan",
    breaks = 10,
    las = 1
    )

### SLIDE 44 ######################################################################################
boxplot(cars)

### SLIDE 45 ######################################################################################
?par
?pdf

### SLIDE 46 ######################################################################################
## Missing values
val <- c(1,3,5,NA,3,6)
val
is.na(val)
which(is.na(val))

### Solution Exercise III - SLIDE 47 ###############################################################
### SLIDE 48 ######################################################################################
val <- c(1,3,5,NA,3,6)
sum(val, na.rm = TRUE)

mean(val)
mean(val, na.rm = TRUE)

### SLIDE 49 ######################################################################################
?save
save(cars.as.list, file="my_cars_as_list.RData")
load(file="my_cars_as_list.RData") #if the file is present in the working directory, if not indicate the path of the .RData file

save(cars.as.list,numeric.vector,rbind.together, file="my_Objects_Rreview_May2016.RData")
load(file="my_Objects_Rreview_May2016.RData") # This will load the 3 objects saved above 

## alternatively
save(cars.as.list, file="my_cars_as_list.RData")
rm(cars.as.list) # To remove the cars.as.list object
ls() # To list all the objects present in the current R session, checking that the cars.as.list is not present anymore
load(file="my_cars_as_list.RData")
ls() # To list all the objects present in the current R session, checking that the cars.as.list is now present

## To save all the objects in the R session
save.image(file="Rreview_2016.RData")
## after closing your R session for example, load the data with:
load(file="Rreview_2016.RData") 

### SLIDE 50 ######################################################################################
## Packages
## install.packages(”PackageName”)
install.packages("heatmap.plus")
## load the library
library("heatmap.plus")
chooseCRANmirror()

### SLIDE 51 ######################################################################################
source("http://bioconductor.org/biocLite.R")
## biocLite(”PackageName”)
biocLite("DESeq2")
## Load the library
library("DESeq2")

#End ######################################################################################


## Extra Slides/code




### Let’s Try Plotting! (IV) - SLIDE 55 ###############################################################
### SLIDE 56 ######################################################################################
# Do a scatter plot with connected dots
plot(x=cars$speed, y=cars$dist,
	xlab = "Speed", ylab = "Distance",
	cex.lab = 1.5,
	main = "A nice scatter plot",
	pch = 16,
	bty = "n",
	col = "dark blue",
	type = "b",
	las = 1
	)

### SLIDE 57 ######################################################################################
# Make your customized version of the boxplot
boxplot(cars,
	width = c(3,1),
	col = "red",
	border = "dark blue",
	names = c("Speed", "Distance"),
	main = "My boxplot",
	notch = TRUE,
	horizontal = TRUE
	)

### SLIDE 58 ######################################################################################
# How can you change the 1.5*IQR parameter?
boxplot(cars,
	width = c(3,1),
	col = "red",
	border = "dark blue",
	names = c("Speed", "Distance"),
	main = "My boxplot",
	range = 1, # this one here!
	notch = TRUE,
	horizontal = TRUE
	)

### SLIDE 59 ######################################################################################
## Print the scatter plot and the boxplot on top of each other and save the figure in a pdf file
pdf("myfigure.pdf", height=10, width=6)
par(mfrow=c(2,1))
plot(x=cars$speed, y=cars$dist,
	xlab = "Speed", ylab = "Distance",
	cex.lab = 1.5,
	main = "A nice scatter plot",
	pch = 16,
	bty = "n",
	col = "dark blue",
	type = "b",
	las = 1
	)

boxplot(cars,
	width = c(3,1),
	col = "red",
	border = "dark blue",
	names = c("Speed", "Distance"),
	main = "My boxplot",
	notch = TRUE,
	horizontal = TRUE
	)
dev.off()


### SLIDE 60 ######################################################################################
x <- 2
if (x>0) { 
  cat("Positive value:",x,"\n")
} else if (x<0) {
  cat("Negative value:",x,"\n")
}

x <- -3
if (x>0) { 
  cat("Positive value:",x,"\n")
} else if (x==0) {
  cat("Zero:",x,"\n")
} else {
  cat("Negative value:",x,"\n")
}


### SLIDE 61 ######################################################################################
for (i in 1:5) {
  cat(i)
}

values <- c(2,1,3,6,5)
for (value in values) {
  ## cat(value)
  print(value)
} 


### SLIDE 62 ######################################################################################
x <- 4
while (x>0) {
 cat("positive value:",x,"\n")
 x <- x-1
}
## Infinite loop!!!
## x <- 4
## while (x>0) {
##  cat("positive value:",x,"\n")
##  x <- x+1
## } 

### Solution exercises V - SLIDE 63 ##############################################################
### SLIDE 64 ######################################################################################

# Print all numbers from 1 to 10
sample <- 1:10 # create the sample to browse

for (n in sample) {
	print (n) # print each number in the sample
	}

# Print all even numbers from 1 to 10
sample <- 1:10 # create the sample to browse

for (n in sample) {
	if (n %% 2 == 0) { # test the rest of the division by 2 (see if even)
		print (n)
		} # no need for a else here (and it is not required)
	}

### SLIDE 65 ######################################################################################
# Print the speed of the first 8 cars using while
# we initialize the index to track how many cars were printed
# 1 to start at the first car
index <- 1								

# we continue until the index is 8
while (index <= 8) {								
	# we access to the speed of the car
	speed <- cars[index, "speed"]
	
	# we print with \n to go to the next line
	cat("Car #", index, "speed:", speed, "\n")
	
	# every iteration we go to the next car
	index <- index + 1
	}

### SLIDE 66 ######################################################################################
# Print the first 8 cars that have a speed more than 9
# we initialize the variables so that they can be used in the loop
counter <- 0 # 0 on the counter to add up
# each time we find an appropriate car, 1 on the index to start at the first car
index   <- 1

while (counter < 8) {
	# to access the speed of the car
	speed <- cars[index, "speed"]
	# if it is more than 9
	if (speed > 9) {
		# we print the car found, with \n to go to the next line 
		cat("Car #", index, "speed:", speed, "\n")
		# and track that we have printed one more car
		counter <- counter + 1
		}
	# every time we go to the next car
	index <- index + 1
	}

### SLIDE 67 ######################################################################################
# How many cars have a speed greater than 10?
# What is their mean distance?
# we put 0 to add up the values we have while browsing the cars
count    <- 0
distance <- 0

# we browse all cars by their index
for (i in 1:nrow(cars)) {
	# test if the speed of the car exceeds 10
	if (cars[i, "speed"] > 10) {
		# here we add it to the group of cars considered
		count    <- count + 1
		# we add its distance to compute the mean afterwards
		distance <- distance + cars[i, "dist"]
		}
	}

# we compute the distance in the end with the global sum of all distances
# and the number of cars used to get that global sum
distanceMean <- distance / count

# print the results, \n is used to go to the next line
cat(count, "cars, mean distance", distanceMean, "\n")

### SLIDE 68 ######################################################################################
function.example <- function(vector.of.values){
  sum.exponent.value <- sum(vector.of.values)^2
  return(sum.exponent.value)
}

dataset.a
function.example(dataset.a)

### SLIDE 69 ######################################################################################
function.example <- function(vector.of.values, exponent.value = 2){
  sum.exponent.value <- sum(vector.of.values)^exponent.value
  return(sum.exponent.value)
}
dataset.a
function.example(dataset.a)
function.example(dataset.a, exponent.value = 10)

### Solution exercises VI - SLIDE 70 ###############################################################
### SLIDE 71 ######################################################################################
# Create a function that takes in a vector and returns its mean.
calculate.mean <- function(x){
	to.return <- mean(x)
	return(to.return)
	}

dataset.a <- c(1, 22, 3, 4, 5)
dataset.a

calculate.mean(dataset.a)

### SLIDE 72 ######################################################################################
# Create a function that takes in a numeric vector and minimum cutoff value. Return the mean,
# median and variance for the numbers in the vector that are greater than the minimum cutoff. Use
# all positive values if the user does not input a minimum cutoff value.
summary.selection <- function(vector.of.values, cutoff.value = 0){
	# select values greater than the cut-off value
	selected     <- vector.of.values[vector.of.values > cutoff.value]
	
	# compute the mean, median, and variance of the selected values
	mean.value   <- mean(selected)
	median.value <- median(selected)
	var.value    <- var(selected)
	
	# return the mean, median, and variance
	to.return    <- list( mean = mean.value, median = median.value, var = var.value)
	return(to.return)
	}

summary.selection(dataset.a)





### END ######################################################################################
##############################################################################################
