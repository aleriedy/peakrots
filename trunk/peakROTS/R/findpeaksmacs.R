
findpeaksmacs.dispatch <- function(arglist) {

	sample <- arglist[4]
	is.original <- (sample == "ORIG") 
	findpeaksmacs(treatment.file = arglist[1], control.file = arglist[2], outputName = arglist[3], sample = sample, is.original)
}

findpeaksmacs <- function(treatment.file, control.file, outputName, sample="", is.original) {

	dir.create(file.path(path.peaks, "log"), showWarnings=FALSE, recursive=TRUE)

	force_max_dup_tags <- !is.original # force max dups for bootstrap/randomized data


        for(shiftsize.index in 1:length(shiftsize)) { # shiftSize
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


## addArgs is for switches (version, help, wig, nolambda, nomodel, diag), e.g.,
## addParams == "--help" returns macs help.
## Parameter names in the params list must be correct parameter names of the MACS.
## See macs --help for complete list of possible parameters.
## stderr of the MACS is redirected to logFile. Use tail -f logFile to follow the output on real time

## macs -- Model-based Analysis for ChIP-Sequencing
##
## Options:
##   -t TFILE, --treatment=TFILE
##                         ChIP-seq treatment files. REQUIRED. When ELANDMULTIPET
##                         is selected, you must provide two files separated by
##                         comma, e.g.
##                         s_1_1_eland_multi.txt,s_1_2_eland_multi.txt
##   -c CFILE, --control=CFILE
##                         Control files. When ELANDMULTIPET is selected, you
##                         must provide two files separated by comma, e.g.
##                         s_2_1_eland_multi.txt,s_2_2_eland_multi.txt
##   --name=NAME           Experiment name, which will be used to generate output
##                         file names. DEFAULT: "NA"
##   --format=FORMAT       Format of tag file, "ELAND" or "BED" or "ELANDMULTI"
##                         or "ELANDMULTIPET". The only acceptable ELAND format
##                         is defined in 00README file, please check it before
##                         you choose ELAND/ELANDMULTI/ELANDMULTIPET. DEFAULT:
##                         "BED"
##   --pvalue=PVALUE       Pvalue cutoff for peak detection. DEFAULT: 1e-5
##
## ------------------------
##
##   --gsize=GSIZE         Effective genome size, default:2.7e+9
##   --tsize=TSIZE         Tag size. DEFAULT: 25
##   --bw=BW               Band width. This value is used while building the
##                         shifting model. If --nomodel is set, 2 time of this
##                         value will be used as a scanwindow width. DEFAULT: 300
##   --mfold=MFOLD         Select the regions with MFOLD high-confidence
##                         enrichment ratio against background to build model.
##                         DEFAULT:32
##
##   --nolambda            If True, MACS will use fixed background lambda as
##                         local lambda for every peak region. Normally, MACS
##                         calculates a dynamic local lambda to reflect the local
##                         bias due to potential chromatin structure.
##   --lambdaset=LAMBDASET
##                         Three levels of nearby region in basepairs to
##                         calculate dynamic lambda, DEFAULT: "1000,5000,10000"
##
##
##   --nomodel             Whether or not to build the shifting model. If True,
##                         MACS will not build model. by default it means
##                         shifting size = 100, try to set shiftsize to change
##                         it. DEFAULT: False
##   --shiftsize=SHIFTSIZE
##                         The arbitrary shift size in bp. When nomodel is true,
##                         MACS will regard this value as 'modeled' d. DEFAULT:
##                         100
##
##
## ------------------------
##   --version             show program's version number and exit
##   -h, --help            show this help message and exit.
##   --verbose=VERBOSE     Set verbose level. 0: only show critical message, 1:
##                         show additional warning message, 2: show process
##                         information, 3: show debug messages. DEFAULT:2
##   --wig                 Whether or not to save shifted raw tag count at every
##                         bp into a wiggle file. WARNING: this process is
##                         time/space consuming!!
##   --wigextend=WIGEXTEND
##                         If set as an integer, when MACS saves wiggle files, it
##                         will extend tag from its middle point to a wigextend
##                         size fragment. By default it is modeled d. Use this
##                         option if you want to increase the resolution in
##                         wiggle file. It doesn't affect peak calling.
##   --space=SPACE         The resoluation for saving wiggle files, by default,
##                         MACS will save the raw tag count every 10 bps. Usable
##                         only with '--wig' option.
##   --diag                Whether or not to produce a diagnosis report. It's up
##                         to 9X time consuming. Please check 00README file for
##                         detail. DEFAULT: False
##   --petdist=PETDIST     Best distance between Pair-End Tags. Only available
##                         when format is 'ELANDMULTIPET'. DEFAULT: 200

## New and hopely better version.
## <...> is for Macs parameters. treatment is the only required parameter.
## How to use:
## runMACS(treatment="Illumina_ChIP-Seq_Demo_Data_Johnson_Science_2007/chip1862_hg18_uniq_with_line_numbers.txt", 
##         control="Illumina_ChIP-Seq_Demo_Data_Johnson_Science_2007/mock1862_hg18_uniq_with_line_numbers.txt", 
##         name="MACS_results/Illumina_ChIP-Seq_Demo_Data", format = "ELAND", 
##         verbose=3, logFile="MACS_results/log/Illumina_ChIP-Seq_Demo_Data.log", 
##         nomodel=TRUE, help=F, version=FALSE)
## Use TRUE or T (FALSE or F) for binary parameters (switches)
## If you don't want to set some parameter (i.e. Macs should use the default values)
## set the parameter value to NULL, e.g, pvalue=NULL or just leave it out. Setting the
## value to null is really handy in loops since you don't have to edit the actual
## function call.


runMACS <- function(..., logFile="/dev/null") {
 
	# Parameter values (character vector)
	params <- c(...) ## Nice :)
	if (!is.null(params)) {

		# Flags
		flags <- paste("--", names(params), sep="")

		# Remove the parameter names
		names(params) <- NULL
    
		# Switches that are on (value is "TRUE" or "T")
		switchOnParams <- NULL
		switchOn <- which(params == "T" | params == "TRUE")
		if (!identical(switchOn, integer(0))) {
			switchOnParams <- flags[switchOn]
			params <- params[-switchOn]
			flags <- flags[-switchOn]
		}

		# Switches that are off are ignored (value is "FALSE" or "F")        
		switchOff <- which(params == "F" | params == "FALSE")
		if (!identical(switchOff, integer(0))) {
			params <- params[-switchOff]
			flags <- flags[-switchOff]
		}
    
		# Command
		command <- paste("macs", paste(flags, params, collapse=" "))
		if (!is.null(switchOnParams)) {
			switchOnParams <-  paste(switchOnParams, collapse=" ")
			command <- paste(command, switchOnParams)
		}
    
		# Run macs. Macs writes its output to stderr (stream number 2)
		# &> redirects both stderr and stdout
		# Iterates through mfold values to find low enough that works
		for (mfold in list(32, 24, 16, 8)) {
			system.output <- system(paste(environment.initialiser, command, paste("--mfold=", mfold, sep=""), "2>", logFile))
			if (system.output == 0) {
				break # was succesfull, don't lower mfold value any more
			}
		}
		
		return(invisible(system.output))
	}
}

# Parses shift size that was used to produce the result file
parseMACSParameters <- function(file){
	output<- read.table(file=file, skip=12, header=FALSE, stringsAsFactors=FALSE, comment.char = "A", fill=TRUE)
	output<- as.numeric(output[1,4])
	output<- ceiling(output/2)
	return(output)
}




