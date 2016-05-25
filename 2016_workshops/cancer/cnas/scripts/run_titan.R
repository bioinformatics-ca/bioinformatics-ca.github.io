#!/usr/bin/env Rscript
# Description: Run TITAN on the HCC1395 Exome
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
library("TitanCNA")

# --- Load the parameters ---
numClusters <- 1
params <- loadDefaultParameters(copyNumber = 5, numberClonalClusters = numClusters, symmetric = TRUE)

# --- Load the tumour allele counts ---
tumAlleleCounts <- loadAlleleCounts("titan/bcftools/tables/HCC1395_exome_tumour.var.het.table.txt")

# For whole exome sequencing data, there will be regions of the genome with no coverage. 
# To account for this, users can indicate the region boundaries (e.g. exons) for which read coverage is expected.
exomeCaptureSpaceDf <- read.table("ref_data/NimbleGenExome_v3.bed", sep = "\t", quote = "", as.is = TRUE)
exomeCaptureSpaceDf <- exomeCaptureSpaceDf[, 1:3]
exomeCaptureSpaceDf[, 1] <- gsub("chr", "", exomeCaptureSpaceDf[, 1])

# --- Correct GC content and mappability biases ---
gcWig <- "ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.gc.wig"
mapWig <- "ref_data/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.map.ws_1000.wig"
tumWig <- "hmmCopy/wig/HCC1395_exome_tumour.wig"
normWig <- "hmmCopy/wig/HCC1395_exome_normal.wig"
cnData <- correctReadDepth(tumWig, normWig, gcWig, mapWig, genomeStyle = "NCBI", targetedSequence = exomeCaptureSpaceDf)

# -- Assign copy number to each position ---
logR <- getPositionOverlap(tumAlleleCounts$chr, tumAlleleCounts$posn, cnData)
tumAlleleCounts$logR <- log(2^logR)
rm(logR, cnData)

# Filter the data
tumAlleleCounts <- filterData(tumAlleleCounts, c(1:22,"X"), minDepth=10, maxDepth=200)

convergeParams <- runEMclonalCN(tumAlleleCounts, gParams = params$genotypeParams, nParams = params$normalParams,
                                pParams = params$ploidyParams, sParams = params$cellPrevParams,
                                maxiter = 20, maxiterUpdate = 1500, txnExpLen = 1e15, txnZstrength = 1e5, 
                                useOutlierState = FALSE, 
                                normalEstimateMethod = "map", estimateS = TRUE, estimatePloidy = TRUE)

optimalPath <- viterbiClonalCN(tumAlleleCounts, convergeParams)

# --- Print out results ---
resultsDir <- "results/titan"
resultFile <- paste(resultsDir, "/HCC1395_exome_tumour.results.txt", sep = "")
dir.create(resultsDir, recursive = TRUE) 
results <- outputTitanResults(tumAlleleCounts, convergeParams, optimalPath, 
                              filename = resultFile, posteriorProbs = FALSE, subcloneProfiles = TRUE)

# -- Print out final model parameters and summary --
outParamFile <- paste(resultsDir, "/HCC1395_exome_tumour.params.txt", sep = "")
outputModelParameters(convergeParams, results, outParamFile)

# --- Plot the LRR and BAF plots for each chromosome ---
ploidy <- convergeParams$phi[length(convergeParams$phi)]

for (chr in c(1:22)) {
  outPlot <- paste(resultsDir, "/", "HCC1395_exome_tumour.", chr, ".png", sep = "")
  png(outPlot, width = 1200, height = 1000, res = 100)
  par(mfrow = c(2,1))
  plotCNlogRByChr(results, chr, ploidy = ploidy, geneAnnot = NULL, spacing = 4, ylim = c(-4,6), cex = 0.5, main = paste("Chromosome", chr, "- LRR"))
  plotAllelicRatio(results, chr, geneAnnot = NULL, spacing = 4, ylim = c(0,1), cex = 0.5, main = paste("Chromosome", chr, "- Allelic Ratio"))
  dev.off()
}
