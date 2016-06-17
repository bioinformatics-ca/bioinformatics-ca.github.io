---
layout: post2
permalink: /analysis_of_metagenomic_data_integrated_assignment1_2016/
title: Analysis of Metagenomic Data 2016 Student Page
header1: Analysis of Metagenomic Data 2016
header2: Integrated Assignment 1
image: CBW_Metagenome_icon.jpg
---

Integrated Assignment 1 Answer Key
----------------------------------

1) The default p-value is 0.01. This can be determined by looking at the help dialogue of the PEAR command, `pear -h`.

2) Since no "step2_otus" or "step3_outs" folders were created, it doesn't look like those steps were performed. We can look at exactly what was done using the log_[...].txt file. In this case, only step 1 and step 4 were performed, which means that we: (A) assigned reads to existing OTUs (closed-reference), (B) took *all* of the unassigned reads and clustered them de novo. This probably took longer than the strategy described in the documents, but it's a good demonstration that you need to keep on top of what an automated pipeline is doing and make sure it's performing the steps you expect.  

3) Check out the files: 
clustering/step1_otus/combined_seqs_otus.log
clustering/step1_otus/sortmerna_otus.log 
- Left over for de novo clustering: 39,788

3) Datasets: SILVA, CORE Oral, Human Oral Microbiome Database (HOMD), NCBI's 16S dataset, or PhytoREF for plastidial rRNA genes. There are many additional datasets for other marker genes (18S, ITS, mitochondrial genes, etc.) Other classification methods include BLAST, rtax, utax, and "16S Classifier". Some reference datasets may be better for specialized environments (such as CORE or HOMD for human oral samples). Every dataset has different requirements before a sequence is included, so be sure to investigate how the reads in your dataset were annotated. Also check the release date for the reference dataset you want to use. GreenGenes, for example, has not been updated since August 2013. You may wish to choose a more frequently updated dataset, like SILVA (last updated April 2015). If you have an organism or group of interest, it is important to check that the organism is well represented in your reference set.

4) Looking at the log script: clustering/log_[...].txt we can see that a script called align_seqs.py was used. Using the documentation for the script: [](http://qiime.org/scripts/align_seqs.html) we can see that the method used is called PyNAST, which is a template-based method. 
FastTree was used to create the phylogenetic tree (see [](http://qiime.org/scripts/make_phylogeny.html)). With the command `FastTree -expert`, you can see that the version used is FastTree 2.1.3 SSE3. This version was released in April 2010. The current version is 2.1.9, released March 2016. Since it often takes developers of pipelines a significant amount of time to update to the most recent version of software tools, it is sometimes prudent to learn to use a tool independent of the pipeline software. This ensures that your results are not needlessly affected by bugs that have already been fixed.

5) The tab-separated format is human readable. In particular, it can be opened and explored in Excel (if the file is small enough) or with programs like `awk` or `grep`. Some programs may not be compatible with the BIOM format, but may be compatible with the simpler format instead.

6) Rarefaction will "throw out" data from samples that have a higher sampling depth than the minimum. If some samples are sequenced to a significantly higher depth, there can be a lot of variation in these samples between independent rarefactions. Often a "jackknifing" procedure is done where the rarefaction is performed many times, allowing a confidence interval to be produced for each sample. Alternatives include variance stabilizing transformations, such as those available in the R DESeq2 and Metagenomeseq packages. A QIIME wrapper script, normalize\_table.py, makes using these packages easier (see <http://qiime.org/scripts/normalize_table.html>).

7) Approximately 30% of the variance was captured by the primary axis, and 13% by the secondary axis, for a total of 45% on a 2D plot. Latitude may be a significant factor affecting the Bray-Curtis distance between points. The ARCT (Arctic) samples mostly seem to group together, but some NWCS (North West Continental Shelf) samples are similar to the Arctic samples.

8) We can open the tab-separated (.tab) OTU tables and search them with a text editor or Excel. Alternatively, we could run `grep Sediminicola otu_table/otu_table.tab` to get the OTU number (the value in the first column). Next, we can grab the sequence from the representative sequence file (rep\_set.fa) with `grep -A 4 '>x$' cluster/rep_set.fasta` where x is replaced by the OTU number found previously (102274), the `-A 4` flag gives us the 4 lines after the matched line. Aligning this sequence against NCBI's database with BLAST shows hits to *Flavobacteria* at 99% identity, but not *Sediminicola*. Both of these genera belong to the Flavobacteriaceae family. It is possible that these unclassified *Flavobacteria* sequences are *Sediminicola*, but it would be healthy to be suspicious of this classification.
