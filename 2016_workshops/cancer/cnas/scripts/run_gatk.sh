#!/bin/bash
# This file will analyse chromosome 21 of the normal and tumour BAM files.
# It will output any position with a variant genotype in one of the two samples.
# The output file can be processed to build the APOLLOH input.
# The output file can also be used to do naive somatic mutation calling.
# WE WILL NOT RUN THIS DURING THE TUTORIAL BECAUSE IT IS A BIT SLOW

ref=/home/ubuntu/CourseData/CG_data/Module5/genome/Homo_sapiens_assembly19.fasta # Reference genome

normalBam=/home/ubuntu/CourseData/CG_data/TCGA/HCC1143/G15511.HCC1143_BL.1.chr21.bam # Normal BAM File

tumourBam=/home/ubuntu/CourseData/CG_data/TCGA/HCC1143/G15511.HCC1143.1.chr21.bam # Tumour BAM File

outFile=/home/ubuntu/CourseData/CG_data/Module5/data/vcf/HCC1143.GATK.vcf 

######### EXECUTE GATK UNIFIED GENOTYPER ###################################################
######### Requires Java to be installed ####################################################
######### -R specifies the reference genome fasta file #####################################
######### -T specifies the tool we want to execute; we want UnifiedGenotyper ###############
######### -L specifies the region we want to analyse. We just pass 21 to specify all of chromosome 21
######### -I specifies input BAM filename
######### -baq RECALCULATE specifies we want to recalculate the base scores
java -jar /usr/local/GATK/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -baq RECALCULATE -L 21 -I $normalBam -I $tumourBam -o $outFile
