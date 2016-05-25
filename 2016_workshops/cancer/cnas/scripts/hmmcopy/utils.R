# Functions for manipulating the output of HMMcopy.
# 
# Author: Andrew Roth
###############################################################################

# Prevent scientific notation in output files. This can cause problems for some 
# software and when displaying in spreadsheets
options(scipen=0)

createResultsDataFrame <- function(sample.id, read.count.data, segments){
	# Converts the results from HMMcopy to data frame which can be written to
	# file. The data frame will have one row per bin.
	
	# Human readable names for the states in the HMM
	state.names <- c('HOMD', 'HETD', 'NEUT', 'GAIN', 'AMP', 'HLAMP')
	
	# Copy state probabilities to a better name and transpose to fit data frame
	state.probs <- t(segments$rho)
	
	# Extract various parameters and format them into a data frame
	results <- cbind(sample = as.character(sample.id), 
			chr = as.character(space(read.count.data)),
			start = start(read.count.data), 
			end = end(read.count.data), 
			state = segments$state,
			event = state.names[segments$state], 
			copy = round(read.count.data$copy, digits = 4),
			pHOMD = sprintf('%0.6f', state.probs[,1]),
			pHETD = sprintf('%0.6f', state.probs[,2]),
			pNEUT = sprintf('%0.6f', state.probs[,3]),
			pGAIN = sprintf('%0.6f', state.probs[,4]),
			pAMP = sprintf('%0.6f', state.probs[,5]),
			pHLAMP = sprintf('%0.6f', state.probs[,6]))
	
	return(results)
}

createSegmentsDataFrame <- function(sample.id, read.count.data, segments){
	# Converts the results from HMMcopy to data frame which can be written to
	# file. The data frame will have one row per segment.
	
	# Round the median value of the segments coverage to 6 digits.
	segments$segs$median <- round(segments$segs$median, digits = 6)
	
	# Compute the bin size
	bin.size <- end(read.count.data[1, ]) - start(read.count.data[1, ]) + 1
	
	# Compute how many bins are in the segments.
	markers <- (segments$segs$end - segments$segs$start + 1) / bin.size
	
	# Human readable names for the states in the HMM
	state.names <- c('HOMD', 'HETD', 'NEUT', 'GAIN', 'AMP', 'HLAMP')
	
	# Convert the column of numeric state values to human readable form
	events <- state.names[as.numeric(as.character(segments$segs$state))]
	
	# Create the data frame
	results <- cbind(sample = as.character(sample.id), 
			segments$segs[, 1:4],
			event = events, 
			bins = markers,
			median = segments$segs$median)
	
	return(results)
}

writeBinLevelResults <- function(sample.id, read.count.data, segments, out.file){
	# Writes the bin level output of HMMcopy to a file. Each row in the output
	# file will correspond to one bin in HMMcopy.
	#
	# Args:
	#	sample.id : Name of sample the data came from
	#	read.count.data : RangedData object with corrected counts.
	#	segments : Segments predicted by HMMcopy
	#	out.file : Path where output file will be written.
	
	bin.results.table <- createResultsDataFrame(sample.id, 
			read.count.data, 
			segments)
	
	write.table(bin.results.table, 
			file=out.file, 
			quote=FALSE, 
			sep="\t",
			row.names=FALSE)
}

writeSegmentResults <- function(sample.id, read.count.data, segments, out.file, 
		igv.compatible=TRUE){
	# Writes the bin level output of HMMcopy to a file. Each row in the output
	# file will correspond to one bin in HMMcopy.
	#
	# Args:
	#	sample.id : Name of sample the data came from
	#	read.count.data : RangedData object with corrected counts.
	#	segments : Segments predicted by HMMcopy
	#	out.file : Path where output file will be written.
	#	igv.compatible : If TRUE will write output in format compatible with IGV
	#					 If FALSE will use format for APOLLOH.
	
	# Get the segment results
	segment.results.table <- createSegmentsDataFrame(sample.id, 
			read.count.data, 
			segments)
	
	if(!igv.compatible){
		# For APOLLOH compatibility we need to re-order the columns 
		segment.results.table <- segment.results.table[, c(1, 2, 3, 4, 7, 8, 5, 6)]
		
		# We also rename the columns
		colnames(segment.results.table) <- c('ID', 
				'chrom', 
				'start', 
				'end', 
				'num.mark', 
				'seg.mean', 
				'state.num', 
				'state.name')
	}
	
	# Output segments to tab separated file
	write.table(segment.results.table, 
			file = out.file,
			quote = FALSE, 
			sep = '\t',
			row.names = FALSE)
}

