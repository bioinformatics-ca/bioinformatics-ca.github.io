# Helpful scripts for running commands and working with the file system.
# 
# Author: Andrew Roth
###############################################################################

safe_make_dir <- function(dir){
	# Creates a directory (including parent directories) if it does not exist.
	if (!file.exists(dir))
	{
		dir.create(dir, recursive=TRUE)
	}
}


