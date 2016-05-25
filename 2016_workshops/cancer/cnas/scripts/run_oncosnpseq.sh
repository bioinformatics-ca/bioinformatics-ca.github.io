#!/bin/bash

TUMOUR_FILE=$1

NORMAL_FILE=$2

SAMPLE_NAME=$3

OUTPUT_DIR=$4

#### OncoSNP Parameter Files ####

STATES_FILE=$ONCOSNPSEQ_DIR/config/tumourstates.txt

HGTABLES_FILE=$ONCOSNPSEQ_DIR/config/hgTables_b37.txt

############## EXECUTE ONCOSNP ####################################################################################
############## For more information, see https://sites.google.com/site/oncosnpseq/usage                 ###########
###################################################################################################################

$ONCOSNPSEQ_DIR/executables/run_oncoseq.sh $MCR_DIR \
    --tumourstatesfile $STATES_FILE \
    --hgtable $HGTABLES_FILE \
    --gcdir $GC_DIR \
    --infile $TUMOUR_FILE \
    --normalfile $NORMAL_FILE \
    --samplename $SAMPLE_NAME \
    --seqtype 'illumina' \
    --outdir $OUTPUT_DIR \
    >$OUTPUT_DIR/run.log 2>$OUTPUT_DIR/run.err
