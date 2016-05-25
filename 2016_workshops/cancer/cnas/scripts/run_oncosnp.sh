#!/bin/bash

BATCH_FILE=$1 

OUTPUT_DIR=$2

#### OncoSNP Parameter Files ####

HYPERPARAMS_FILE=$ONCOSNP_DIR/data/cyau/temp/oncosnp/configuration/hyperparameters-affy.dat

LEVELS_FILE=$ONCOSNP_DIR/configuration/levels-affy.dat

TRAINING_FILE=$ONCOSNP_DIR/configuration/trainingStates.dat

STATES_FILE=$ONCOSNP_DIR/configuration/tumourStates.dat

HGTABLES_FILE=$ONCOSNP_DIR/configuration/hgTables_b37.txt

############## EXECUTE ONCOSNP ####################################################################################
############## For more information, see https://sites.google.com/site/oncosnp/user-guide/using-oncosnp ###########
###################################################################################################################

#Arguments:
#--batch-file {filename}
#Specifies the location of the batch file containing the list of files to process.

#--output-dir {directory}
#Specifies the directory where output files will be placed.

#--gcdir {directory}
#Directory where the local GC content files are contained. Specifying this option activates the local GC content based correction which corrects for the "wave effect" often observed in the Log R Ratio.

#--paramsfile {filename}
#Location of the parameters file.

#--levelsfile {filename}
#Location of the configuration file containing Log R Ratio levels for different copy number states. 

#--trainingstatesfile {filename}
#Location of configuration file containing short list of tumour states used for training. 

#--tumourstatesfile {filename}
#Location of configuration file containing full list of tumour states.

#--subsample {number}
#Uses a subset of the data for training and learning the noise model. A value of 10 would use 1/10th of the data for training. We recommend a suitable value that gives a training dataset of between 20-30,000 probes for reasonable computational times. 

#--emiters {number}
#Maximum number of expectation-maximisation steps to use during training.

#--stromal
#Activates the stromal contamination model. The program will attempt to estimate the level of stromal contamination in the sample (up to 50%).

#--intratumour
#Activates the intra-tumour heterogeneity model. The program will attempt to estimate the level of intra-tumour heterogeneity in the sample.

#NOTE: When both stromal contamination and intra-tumour heterogeneity models are activated, the program performs a joint analysis and attempts to estimate the baseline level of stromal contamination and any intra-tumour heterogeneity occurring in addition to this. This is the most computationally demanding mode of operation but possibly also the most accurate.

#--chr {number}
#Processed a specified set of chromosomes only, e.g. --chr 1,2,5:12, would analyse chromosomes 1, 2 and 5 to 12 only. NOTE: Currently Chromosome MT probes are not processed.

#--female
#Adds the X chromosome to the list of chromosomes to be processed.

#--fulloutput
#Generates a probe-by-probe classification summary file.

#--quickmode 
#Use a quicker estimation (possibly less accurate) method for fitting the OncoSNP model to the sample data.

$ONCOSNP_DIR/run_oncosnp.sh $MCR_DIR \
    --batch-file $BATCH_FILE \
    --output-dir $OUTPUT_DIR \
    --fulloutput --plot \
    --gcdir $GC_DIR \
    --paramsfile $HYPERPARAMS_FILE \
    --levelsfile $LEVELS_FILE \
    --subsample 30 \
    --emiters 1 \
    --female \
    --trainingstatesfile $TRAINING_FILE \
    --tumourstatesfile $STATES_FILE \
    --chr 21 \
    --stromal \
    --intratumor \
    --hgtables $HGTABLES_FILE \
    >$OUTPUT_DIR/run.log 2>$OUTPUT_DIR/run.err





