---
layout: post2
permalink: /BiCG_Module8_Annovar_lab/
title: BiCG Module 8 Annovar Lab
header1: BiCG Module 8 Annovar Lab
header2: Workshop pages for students
image: CBW_cancerDNA_icon-16.jpg
---

## Setup

First login into the server.

~~~bash
ssh -i CBW.pem ubuntu@cbw##.dyndns.info
~~~

Now enter the ~/workspace directory

~~~bash
cd ~/workspace
~~~

Create a directory for this module and enter this directory:

~~~bash
mkdir Module8
cd Module8
~~~

Download the Data Set Input - VCF file to your local computer.  In a separate local machine terminal instance, move into the directory where you have downloaded the VCF to and copy the VCF to the workspace.

~~~bash
scp -i CBWCG.pem ./passed.somatic.snvs.vcf ubuntu@cbw##.dyndns.info://home/ubuntu/workspace/Module8
~~~

Convert the VCF file to an annovar input file.

~~~bash
convert2annovar.pl --includeinfo -format vcf4old Passed.somatic.snvs.vcf > Passed.somatic.snvs.vcf.annovar.in.txt
~~~

The output you will see is:

~~~bash
NOTICE: Read 121 lines and wrote 0 different variants at 8 genomic positions (8 SNPs and 0 indels)
NOTICE: Among 8 different variants at 8 positions, 0 are heterozygotes, 0 are homozygotes
NOTICE: Among 8 SNPs, 1 are transitions, 7 are transversions (ratio=0.14)
~~~

Run table_annovar.pl to annotate the variants in the annovar input file you have created.

~~~bash
table_annovar.pl --buildver hg19 Passed.somatic.snvs.vcf.annovar.in.txt /media/cbwdata/software/annovar/annovar/humandb/ --protocol refGene,ljb26_all,1000g2014oct_all,caddgt10,cg69,clinvar_20150330,cosmic70,esp6500siv2_all,exac02,snp138,genomicSuperDups,phastConsElements46way --operation g,f,f,f,f,f,f,f,f,f,r,r --nastring NA --outfile passed.somatic.snvs.vcf.annovar.out.txt
~~~

#### Arguments in the command


**--buildver hg19**   
* human genome reference build

**Passed.somatic.snvs.vcf.annovar.in.txt** 
* input file

**/media/cbwdata/software/annovar/annovar/humandb/** 
* Annovar db path

**--protocol refGene,ljb26_all,1000g2014oct_all,caddgt10,cg69,clinvar_20150330,cosmic70,esp6500siv2_all,exac02,snp138,genomicSuperDups,phastConsElements46way**
* list of Annovar annotation modules to be executed, corresponding to specific databases

**--operation g,f,f,f,f,f,f,f,f,f,r,r**
* type of operation to be executed by Annovar annotation modules
   * g = gene (only for the gene database), f = filter (exact match by coordinates, ref, alt), r = regional (coordinate overlap)
* --protocol and --operation need the same number of comma-separated items

**--nastring NA**
* encoding of NA values, "NA" is good for R post-processing, use "." for VCF output

**--outfile passed.somatic.snvs.vcf.annovar.out.txt**
* output file name prefix

#### Output File

The output file created (Passed.somatic.snvs.vcf.annovar.out.txt.hg19_multianno.txt) is a tab-delimited file, where each row represents one variant, and each column represents one annotation task. Table_annovar allows you to specify exactly which columns or annotation tasks are required, and allows you to select multiple versions of the same analysis (such as multiple gene-definition systems or multiple dbSNP databases).

The output you will see is:

~~~bash
 NOTICE: Processing operation=g protocol=refGene

 NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg19 -dbtype refGene -outfile passed.somatic.snvs.vcf.annovar.out.txt.refGene -exonsort passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: Reading gene annotation from /usr/local/annovar/humandb/hg19_refGene.txt ... Done with 51039 transcripts (including 11569 without coding sequence annotation) for 26311 unique genes
 NOTICE: Reading FASTA sequences from /usr/local/annovar/humandb/hg19_refGeneMrna.fa ... Done with 21 sequences
 WARNING: A total of 345 sequences will be ignored due to lack of correct ORF annotation
 NOTICE: Finished gene-based annotation on 8 genetic variants in passed.somatic.snvs.vcf.annovar.in.txt
 NOTICE: Output files were written to passed.somatic.snvs.vcf.annovar.out.txt.refGene.variant_function, passed.somatic.snvs.vcf.annovar.out.txt.refGene.exonic_variant_function
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=ljb26_all
 NOTICE: Finished reading 25 column headers for '-dbtype ljb26_all'

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype ljb26_all -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/ -otherinfo>
 NOTICE: the --dbtype ljb26_all is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_ljb26_all_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_ljb26_all_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 557362 and the number of bins to be scanned is 7
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_ljb26_all.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=1000g2014oct_all

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2014oct_all -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_ALL.sites.2014_10_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_ALL.sites.2014_10_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 2824642 and the number of bins to be scanned is 6
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_ALL.sites.2014_10.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=caddgt10

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype caddgt10 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype caddgt10 is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_caddgt10_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_caddgt10_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 2625942 and the number of bins to be scanned is 6
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_caddgt10.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=cg69

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype cg69 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype cg69 is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_cg69_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_cg69_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 2789339 and the number of bins to be scanned is 6
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_cg69.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=clinvar_20140929

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype clinvar_20140929 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype clinvar_20140929 is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_clinvar_20140929_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_clinvar_20140929_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 44738 and the number of bins to be scanned is 1
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_clinvar_20140929.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=cosmic70

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype cosmic70 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype cosmic70 is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_cosmic70_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_cosmic70_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 232279 and the number of bins to be scanned is 5
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_cosmic70.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=esp6500siv2_all

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype esp6500siv2_all is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_esp6500siv2_all_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_esp6500siv2_all_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 594771 and the number of bins to be scanned is 7
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_esp6500siv2_all.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=exac02

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype exac02 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: the --dbtype exac02 is assumed to be in generic ANNOVAR database format
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_exac02_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_exac02_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 750585 and the number of bins to be scanned is 7
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_exac02.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=f protocol=snp138

 NOTICE: Running system command <annotate_variation.pl -filter -dbtype snp138 -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: Variants matching filtering criteria are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_snp138_dropped, other variants are written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_snp138_filtered
 NOTICE: Processing next batch with 8 unique variants in 8 input lines
 NOTICE: Database index loaded. Total number of bins is 2894320 and the number of bins to be scanned is 6
 NOTICE: Scanning filter database /usr/local/annovar/humandb/hg19_snp138.txt...Done
 -----------------------------------------------------------------
 NOTICE: Processing operation=r protocol=genomicSuperDups

 NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype genomicSuperDups -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: Reading annotation database /usr/local/annovar/humandb/hg19_genomicSuperDups.txt ... Done with 51599 regions
 NOTICE: Finished region-based annotation on 8 genetic variants in passed.somatic.snvs.vcf.annovar.in.txt
 NOTICE: Output file is written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_genomicSuperDups
 -----------------------------------------------------------------
 NOTICE: Processing operation=r protocol=phastConsElements46wayPlacental

 NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype phastConsElements46wayPlacental -buildver hg19 -outfile passed.somatic.snvs.vcf.annovar.out.txt passed.somatic.snvs.vcf.annovar.in.txt /usr/local/annovar/humandb/>
 NOTICE: Reading annotation database /usr/local/annovar/humandb/hg19_phastConsElements46wayPlacental.txt ... Done with 3743478 regions
 NOTICE: Finished region-based annotation on 8 genetic variants in passed.somatic.snvs.vcf.annovar.in.txt
 NOTICE: Output file is written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_phastConsElements46wayPlacental
 -----------------------------------------------------------------
 NOTICE: Multianno output file is written to passed.somatic.snvs.vcf.annovar.out.txt.hg19_multianno.txt
~~~

From a separate local machine terminal instance, copy the output file back to your local machine

~~~bash
scp -i CBW.pem ubuntu@cbw##.dyndns.info://home/ubuntu/workspace/Module8/passed.somatic.snvs.vcf.annovar.out.txt.hg19_multianno.txt ./
~~~

You can open the file in Excel (select "tab-delimited" when opening the file). Click the "DATA" tab at the menu bar, then click the big "Filter" button. Then click any one of the headings to filter out variants, essentially by clicking the check boxes. 
