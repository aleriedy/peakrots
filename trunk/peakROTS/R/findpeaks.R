
findpeaks.dispatch <- function(arglist) {

	sample <- arglist[4]
	is.original <- (sample == "ORIG") 
	findPeaks(treatment.file = arglist[1], control.file = arglist[2], outputName = arglist[3], sample = sample, is.original)
}

findPeaks <- function(treatment.file, control.file, outputName, sample="", is.original) {

	dir.create(file.path(path.peaks, "log"), showWarnings=FALSE, recursive=TRUE)

	force_max_dup_tags <- !is.original # force max dups for bootstrap/randomized data


-       for(shiftsize.index in 1:length(shiftsize)) { # shiftSize
		for(tsize.index in 1:length(tsize)) { # tSize
	        	for(bw.index in 1:length(bw)) { # bandWidth
	            		for(nolambda.index in 1:length(nolambda)) { # noLambda
	            			for(mfold.index in 1:length(mfold)) { # mfold
					
					do.model.shiftsize <- is.null(shiftsize[[shiftsize.index]])
					if (is.original && do.model.shiftsize) {
						actual.shiftsize = NULL
						nomodel = FALSE

					} else if (!is.original && do.model.shiftsize) {
                            			actual.shiftsize <- parseMACSParameters(file.path(path.peaks, paste("MacspeaksOriginal", 
                                    			"_SHIFT",shiftsize[[shiftsize.index]], 
			                                "_TAG",tsize[[tsize.index]], 
                        			        "_BW",bw[[bw.index]], 
			                                "_NL",nolambda[[nolambda.index]], 
			                                "_MFOLD",mfold[[mfold.index]], 
							"_ORIG_peaks.xls", 
			                                sep="")))
						nomodel = TRUE

					} else {
						actual.shiftsize = shiftsize[[shiftsize.index]]
						nomodel = TRUE
					}

						
					runMACS(
						treatment = treatment.file,
						control = control.file,
						name = file.path(path.peaks, paste(
							outputName, 
							"_SHIFT", shiftsize[[shiftsize.index]], 
							"_TAG", tsize[[tsize.index]], 
							"_BW", bw[[bw.index]], 
							"_NL", nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]], 
							"_", sample, sep="")),	
						format = format,
						verbose = 3,
						logFile = file.path(path.peaks, "log", paste(paste(
							outputName, 
							"_SHIFT", shiftsize[[shiftsize.index]], 
							"_TAG", tsize[[tsize.index]], 
							"_BW", bw[[bw.index]], 
							"_NL", nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]],
							"_", sample, sep=""), "log", sep=".")),


						## Modifications needed ----------------------------------------------
						shiftsize = actual.shiftsize,
						nomodel = nomodel,
						tsize = tsize[[tsize.index]],
						bw = bw[[bw.index]],
						nolambda = nolambda[[nolambda.index]],
						mfold = mfold[[mfold.index]],
						pvalue = pvalue,
						force_max_dup_tags = force_max_dup_tags,
						max_dup_tags = 100
						## -------------------------------------------------------------------
					)
					} # mfold
				} # noLambda
			} # bandWidth
		} # tSize
	} # shiftSize	
}

findPeaks.PeakSeq <- function(treatment.file, control.file, outputName, sample="", is.original) {

	dir.create(file.path(path.peaks, "log"), showWarnings=FALSE, recursive=TRUE)

	force_max_dup_tags <- !is.original # force max dups for bootstrap/randomized data

        for(READ_LENGTH.index in 1:length(READ_LENGTH)) { # parameter: READ_LENGTH
		for(W_SIZE.index in 1:length(W_SIZE)) { # parameter: W_SIZE
	        	for(W_PER_C.index in 1:length(W_PER_C)) { # parameter: W_PER_C
					
					#do.model.shiftsize <- is.null(shiftsize[[shiftsize.index]])
					#if (is.original && do.model.shiftsize) {
					#	actual.shiftsize = NULL
					#	nomodel = FALSE
					#
					#} else if (!is.original && do.model.shiftsize) {
                            		#	actual.shiftsize <- parseMACSParameters(file.path(path.peaks, paste("MacspeaksOriginal", 
                                    	#		"_SHIFT",shiftsize[[shiftsize.index]], 
			                #                 "_TAG",tsize[[tsize.index]], 
                        		#	        "_BW",bw[[bw.index]], 
			                #                 "_NL",nolambda[[nolambda.index]], 
			                #                "_MFOLD",mfold[[mfold.index]], 
					#		"_ORIG_peaks.xls", 
			                #                sep="")))
					#	nomodel = TRUE
					#
					#} else {
					#	actual.shiftsize = shiftsize[[shiftsize.index]]
					#	nomodel = TRUE
					#}

						
					runMACS(
						treatment = treatment.file,
						control = control.file,
						name = file.path(path.peaks, paste(
							outputName, 
							"_SHIFT", shiftsize[[shiftsize.index]], 
							"_TAG", tsize[[tsize.index]], 
							"_BW", bw[[bw.index]], 
							"_NL", nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]], 
							"_", sample, sep="")),	
						format = format,
						verbose = 3,
						logFile = file.path(path.peaks, "log", paste(paste(
							outputName, 
							"_SHIFT", shiftsize[[shiftsize.index]], 
							"_TAG", tsize[[tsize.index]], 
							"_BW", bw[[bw.index]], 
							"_NL", nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]],
							"_", sample, sep=""), "log", sep=".")),
					)

			} # W_PER_C
		} # W_SIZE
	} # READ_LENGTH
}
