# Plotting functions for HMMcopy.
# 
# Author: Andrew Roth
###############################################################################

plotBiases <- function(corrected.data, plot.file){
	# This function will create the density plots for GC and mappability bias.
	#
	# Args:
	#	corrected.data : RangedData object containing corrected read counts.
	#	plot.file : Path where plot file will be written in PDF format. 
	
	# Begin plot to be saved as PDF
	pdf(plot.file, width=15, height=15)
	
	# Create an area to plot the densities
	par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
	
	# Plot the results
	plotBias(corrected.data, pch = 20, cex = 0.5)
	
	# Close the device and write the file
	dev.off()
}

plotReadCounts <- function(corrected.data, plot.file){
	# This function will create the chromosome plots for uncorrected and 
	# corrected read counts.
	#
	# Args:
	#	corrected.data : RangedData object containing corrected read counts.
	#	plot.file : Path where plot file will be written in PDF format.
	# Begin plot to be saved as PDF
	pdf(plot.file, width=15, height=15)
	
	# Create an area to plot the densities
	par(mar = c(4, 4, 2, 0))
	
	# Plot the results
	plotCorrection(corrected.data, pch = '.')
	
	# Close the device and write the file
	dev.off()
}
