############################
#
#   OICR Cancer Workshop Module 9: Data Integration and Survival Workshop
#   Lab Script
#   
#   Lauren Erdman
#   
############################

############################
#
#   NOTE: THIS CODE REQUIRES THE FOLLOWING PACKAGES TO BE INSTALLED
#   SNFtool
#   RColorBrewer
#   survival
#   rms
#
############################

## Command line package installation
install.packages(pkgs = c("SNFtool","RColorBrewer","survival","rms"))

### Set your own working directory
setwd("Your working directory here")

load("OICR-Survival-Workshop-Data-revised-6-3-2016.RData")

##############################
#**********************
#########
#
#     Similarity Network Fusion
#
#########
#**********************
###############################
library('SNFtool')
## First, set all the parameters:
K = 20;    # number of neighbors, usually (10~30)
alpha = 0.5;    # hyperparameter, usually (0.3~0.8)
T = 20; 	# Number of Iterations, usually (10~20)

### 
#   TRANSPOSE AND STANDARDIZE DATA GENOMIC DATA BEING USED
###

std.kirc.list <- lapply(X = list(methyl = t(methyl.kirc),
                                 mirna = t(mirna.kirc),
                                 mrna = t(mrna.kirc)),
                        standardNormalization)

###
#   GENERATE DISTANCE MATRICES USING EUCLIDEAN DISTANCE
###

dist.kirc.matrices <- lapply(X = std.kirc.list,
                             function(x){dist2(as.matrix(x),as.matrix(x))})

###
#   GENERATE AFFINITY MATRICES
###

affinity.kirc.matrices <- lapply(X = dist.kirc.matrices,
                                 function(x){affinityMatrix(x,K,alpha)})

###
#   CLUSTER INDIVIDUAL DATA TYPES  
###

(n.clusters.estimated <- lapply(X = affinity.kirc.matrices,
                               function(x){estimateNumberOfClustersGivenGraph(x)[[1]]}))

clustered.groups <- sapply(X = seq(1,3),
                           function(x){spectralClustering(affinity = affinity.kirc.matrices[[x]],
                                                          K = n.clusters.estimated[[x]])})


colnames(clustered.groups) <- names(affinity.kirc.matrices)

## Looking at distribution of group assignment
apply(clustered.groups,2,table)
    ## Weird! Looks like there is a group of outliers in miRNA and mRNA

## Using heatmap to look at clusters:


displayClustersWithHeatmap(affinity.kirc.matrices[['methyl']],group = clustered.groups[,'methyl'])
displayClustersWithHeatmap(affinity.kirc.matrices[['mirna']],group = clustered.groups[,'mirna'])
displayClustersWithHeatmap(affinity.kirc.matrices[['mrna']],group = clustered.groups[,'mrna'])
## This color isn't the best, let's revise the heatmap colors using RColorBrewer
###
#   CHANGING HEATMAP COLORS USING RColorBrewer
###

library("RColorBrewer")
rd.gy = colorRampPalette(brewer.pal(n = 11,name = "RdGy"))(50)

?colorRampPalette

display.brewer.all()

displayClustersWithHeatmap(affinity.kirc.matrices[['methyl']],group = clustered.groups[,'methyl'],col = rd.gy)
displayClustersWithHeatmap(affinity.kirc.matrices[['mirna']],group = clustered.groups[,'mirna'],col = rd.gy)
displayClustersWithHeatmap(affinity.kirc.matrices[['mrna']],group = clustered.groups[,'mrna'],col = rd.gy)
  ## those groups of 3 seem to be ruining our signal elsewhere - let's see what removing them gets us

####
#   REMOVING INDIVIDUALS FROM EACH DATASET
####  
# First, let's make a vector of IDs we'd like to keep 


ids.to.keep <- colnames(mirna.kirc)[which(clustered.groups[,'mrna'] != 3)]

sum(clinic.kirc$revised.ids %in% ids.to.keep) ## number of individuals we expect in each dataset

####
#   SUBSETTING OUR DATASETS
####

## We handle the clinic data separately since the structure is different
##    than the genomics data


clinic.kirc.sub <- clinic.kirc[match(ids.to.keep,clinic.kirc$revised.ids),]



genomic.sub.list <- lapply(list(methyl=t(methyl.kirc),
                                mirna = t(mirna.kirc),
                                mrna = t(mrna.kirc)),
                           function(x){x[match(ids.to.keep,rownames(x)),]})

### 
#   DOUBLE CHECKING THAT OUR DATA IS NAMED CORRECTLY AND HAS THE CORRECT DIMENSIONS
###
names(genomic.sub.list)
lapply(genomic.sub.list,dim)
str(genomic.sub.list)

#########
#
#   RE-RUNNING DATA-SPECIFIC AFFINITY MATRIX CLUSTERING WITH SUBSET DATA
#
#########

### 
#   STANDARD NORMALIZE DATA GENOMIC DATA BEING USED
###

std.kirc.list.sub <- lapply(X = genomic.sub.list,
                            standardNormalization)

###
#   GENERATE DISTANCE MATRICES USING EUCLIDEAN DISTANCE
###

dist.kirc.matrices.sub <- lapply(X = std.kirc.list.sub,
                                 function(x){dist2(as.matrix(x),as.matrix(x))})

###
#   GENERATE AFFINITY MATRICES
###

affinity.kirc.matrices.sub <- lapply(X = dist.kirc.matrices.sub,
                                     function(x){affinityMatrix(x,K,alpha)})

###
#   CLUSTER INDIVIDUAL DATA TYPES  
###

n.clusters.estimated.sub <- lapply(X = affinity.kirc.matrices.sub,
                                   function(x){estimateNumberOfClustersGivenGraph(x)[[1]]})
clustered.groups.sub <- sapply(X = seq(1,3),
                               function(x){spectralClustering(affinity = affinity.kirc.matrices.sub[[x]],
                                                              K = n.clusters.estimated.sub[[x]])})
colnames(clustered.groups.sub) <- names(affinity.kirc.matrices.sub)

## Looking at distribution of group assignment  
apply(clustered.groups.sub,2,table)
## No obvious outliers! :) 

## Using heatmap to look at clusters:
displayClustersWithHeatmap(affinity.kirc.matrices.sub[['methyl']],
                           group = clustered.groups.sub[,'methyl'],col = rd.gy)
displayClustersWithHeatmap(affinity.kirc.matrices.sub[['mirna']],
                           group = clustered.groups.sub[,'mirna'],col = rd.gy)
displayClustersWithHeatmap(affinity.kirc.matrices.sub[['mrna']],
                           group = clustered.groups.sub[,'mrna'],col = rd.gy)


## Much nicer clusters! 

###
#   RUN SNF
###

kirc.snf <- SNF(affinity.kirc.matrices.sub,K = K,t = T)

## Naming names and columns
colnames(kirc.snf) <- rownames(kirc.snf) <- rownames(genomic.sub.list[[1]])
str(kirc.snf)

###
#   FIND NUMBER OF CLUSTERS
###

estimateNumberOfClustersGivenGraph(W = kirc.snf,NUMC = 2:5)

###
#   GENERATE GROUP ASSIGNMENTS FROM NUMBER OF CLUSTERS DEFINED ABOVE
###

snf.groups <- spectralClustering(kirc.snf,K = 2)

## LOOK AT GROUP SIZES
table(snf.groups)

## SET UP A DATAFRAME WITH GROUP ASSIGNMENT BY ID 
##    (WE WILL USE THIS IN THE SURVIVAL ANALYSIS)
ids.groups2 <- data.frame(cbind(colnames(kirc.snf),snf.groups))
names(ids.groups2) <- c("id","group")

head(ids.groups2)

###
#   GENERATE HEATMAP OF RESULTING MATRIX WITH CLUSTERING
###
displayClustersWithHeatmap(W = kirc.snf,group = snf.groups,col = rd.gy)

##############################
#**********************
#########
#
#     Survival Analysis
#
#########
#**********************
###############################

###########################
#     GENERATING SURVIVAL OUTCOME
###########################

### getting names columns in clinic data (don't forget - it's subsetted now)
names(clinic.kirc.sub)

## first take a look at the variables we're going to use
# str(clinic.kirc$patient.age_at_initial_pathologic_diagnosis)
head(clinic.kirc.sub$patient.days_to_death)
head(clinic.kirc.sub$patient.days_to_last_known_alive)
head(clinic.kirc.sub$patient.days_to_last_followup)
head(clinic.kirc.sub$patient.vital_status)

clinic.kirc.sub$time.to.event <- clinic.kirc.sub$patient.days_to_last_followup
clinic.kirc.sub$time.to.event[is.na(clinic.kirc.sub$patient.days_to_death) == FALSE] <- 
  clinic.kirc.sub$patient.days_to_death[is.na(clinic.kirc.sub$patient.days_to_death) == FALSE]

clinic.kirc.sub$event <- NA
clinic.kirc.sub$event[clinic.kirc.sub$patient.vital_status == "alive"] <- 0
clinic.kirc.sub$event[clinic.kirc.sub$patient.vital_status == "dead"] <- 1

library("survival")
clinic.kirc.sub$survival.outcome <- Surv(clinic.kirc.sub$time.to.event,
                                         clinic.kirc.sub$event)

#######################
#   SUMMARIZING SURVIVAL WITHOUT COVARIATE
#######################
(kirc.survival.fit <- survfit(clinic.kirc.sub$survival.outcome ~ 1, 
                              conf.type = "log-log"))

#######
#     CREATING BASIC KM CURVE
#######
plot(kirc.survival.fit,col="blue4")

#######
#     COMPARING SURVIVAL ACROSS GROUPS
#######

## merging SNF groups and clinic.data (recall we made 'ids.groups2' when we clustered SNF)
clinic.kirc.snf.group <- merge(x = clinic.kirc.sub, ## clinic dataframe
                               y = ids.groups2, ## SNF group dataframe
                               by.x="revised.ids", ## clinic ID column
                               by.y="id")  ## SNF group ID column

## creating factor variables for cluster and sex
clinic.kirc.snf.group$sex <- factor(clinic.kirc.snf.group$patient.gender,
                                    levels = c("male","female"))
clinic.kirc.snf.group$cluster <- factor(clinic.kirc.snf.group$group,
                                        levels = c(1,2))

###
#   TESTING THE DIFFERENCE IN SURVIVAL TIME USING PETO&PETO MODIVICATION ON THE 
#       GEHAN-WILCOXON TEST, USING THE survdiff FUNCTION
###

## Sex

survdiff(clinic.kirc.snf.group$survival.outcome ~ 
           clinic.kirc.snf.group$sex, rho=1)

## Cluster Assignment

survdiff(clinic.kirc.snf.group$survival.outcome ~ 
           clinic.kirc.snf.group$cluster, rho=1)

# npsurv (non-parametric survival fit) function is a work around/replacement 
#   for survfit since survfit no longer works with survplot which we want to use below

library('rms')

## looking at marginal survival difference by sex

kirc.survival.fit.by.sex <- npsurv(survival.outcome ~ sex,
                                   data = clinic.kirc.snf.group)

## looking at marginal survival difference by SNF generated cluster

kirc.survival.fit.by.snf.group <- npsurv(survival.outcome ~ cluster,
                                         data = clinic.kirc.snf.group)

#### 
#   ADDING CONFIDENCE BOUNDS AND COLORS TO KM CURVES PLOTTED BY GROUP
####  

survplot(fit = kirc.survival.fit.by.sex,col=c('blue','deeppink2'),
         col.fill = sapply(c('blue','deeppink2'),function(x){adjustcolor(x, alpha.f = 0.1)}),
         xlab="Days to Event")

survplot(fit = kirc.survival.fit.by.snf.group,col=c('forestgreen','darkorchid4'),
         col.fill = sapply(c('forestgreen','darkorchid4'),function(x){adjustcolor(x, alpha.f = 0.15)}),
         xlab="Days to Event")


###########
#   Cox proportional hazards analysis
###########

# Want to fit model with survival as an outcome, 
  # analyzing cluster assignment while controling for sex as covariates 

## Efron method
(coxph.fit <- coxph(survival.outcome ~ 
                      cluster + sex + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.kirc.snf.group,method = "efron"))
summary(coxph.fit)

## Exact method
(coxph.fit <- coxph(survival.outcome ~ 
                      cluster + sex + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.kirc.snf.group,method = "exact"))
summary(coxph.fit)

## Breslow method
(coxph.fit <- coxph(survival.outcome ~ 
                      cluster + sex + patient.age_at_initial_pathologic_diagnosis, 
                    data = clinic.kirc.snf.group,method = "breslow"))
summary(coxph.fit)

## Extracting results
str(summary(coxph.fit))

(coxph.coefs <- summary(coxph.fit)$coef)
(coxph.confint <- summary(coxph.fit)$conf.int)

(coxph.results <- cbind(coxph.coefs,coxph.confint))
colnames(coxph.results)

write.csv(coxph.results[,c("coef","exp(coef)","se(coef)","Pr(>|z|)","lower .95","upper .95" )],
          file = "Cox-PH-model-results.csv")

###
#   Using cox.zph to test for covariate-specific 
#     and global proportional hazards as well as
#     plotting scho residuals to check for 
#     non-proportional hazards -- significance implies non-proportionality
###

cox.zph(fit = coxph.fit)  
par(mfrow=c(2,2))
plot(cox.zph(fit = coxph.fit))
par(mfrow=c(1,1))

### 
#   Checking for influential observations (outliers)
###

dfbeta <- residuals(coxph.fit, type = 'dfbeta') ## Dataframe of change in coefficients as each individual removed

par(mfrow=c(2,2))
for(j in 1:3){
  plot(dfbeta[,j],ylab=names(coef(coxph.fit))[j],
       pch=19,col='blue')
  abline(h=0,lty=2)
}
par(mfrow=c(1,1))
##  No terribly influential points


###
#   Checking for linearity in the covariates using plots of 
#   martingale residuals against the individual covariates
#     NOTE: This is not necessary for binary variables
#         so we only check it in our age of initial diagnosis 
#           covariate
###
martingale.resids <- residuals(coxph.fit,type = 'martingale')
seq(1,nrow(clinic.kirc.snf.group))[!(seq(1,nrow(clinic.kirc.snf.group)) %in% as.numeric(names(martingale.resids)))] 
  ## checking that there are no missing residual values using the indices of the martingale residuals
    # (the value of 'integer(0)' being returned tells us we haven't missed anything)

par(mfrow=c(2,1))
## Null plot for residuals: 
plot(y = martingale.resids,
     x = clinic.kirc.snf.group$patient.age_at_initial_pathologic_diagnosis,
     ylab = 'Residuals', xlab = 'Age of Initial Diag', pch= 19, col = 'blue')
abline(h=0,lty=2,col='red')
lines(lowess(x = clinic.kirc.snf.group$patient.age_at_initial_pathologic_diagnosis,
             y = martingale.resids, iter = 0))

## Component-plus-residual plot


b <- coef(coxph.fit)[3]
x <- clinic.kirc.snf.group$patient.age_at_initial_pathologic_diagnosis
plot(x, b*x + martingale.resids, 
     xlab='Age of Initial Diag',
     ylab="component+residual", 
     pch = 19, col = 'blue')
abline(lm(b*x + martingale.resids ~ x), lty=2, col = 'red')
lines(lowess(x, b*x + martingale.resids, iter = 0))

  ## deviation of lowess line from 0-line and fit slope are 
    # extremely small therefore linearity seems to hold

par(mfrow=c(1,1))

############################
#
#       BONUS MATERIAL:  
#         Making predictions using our Cox PH model
# 
############################

# Suppose you now have a new individual you'd like to predict the survival of: 
individual_new <- data.frame(cluster=factor(2),sex="male",patient.age_at_initial_pathologic_diagnosis=21)
      ### Note the way each input variable has to be named the EXACT way it was fit in our cox model of fit
                            ## and the cluster value has to be in the "factor" form (set using the factor function)
  
predict(coxph.fit,individual_new,type="risk") ## this is the risk of your 21 year old group 2 male patient 
                                                  ## relative to the average of your sample
