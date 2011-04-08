## PeakSeq ROTS, T.D. Laajala, tlaajala@cc.hut.fi
## Summer 2010, Turku University, Math Department

findpeakspeakseq.dispatch <- function(arglist) {

	sample <- arglist[5]
	is.original <- (sample == "ORIG") 
	findpeakspeakseq(data.path = arglist[1], treatment.file = arglist[2], control.file = arglist[3], outputName = arglist[4], sample = sample, is.original)
}

findpeakspeakseq <- function(data.path, treatment.file, control.file, outputName, sample="", is.original) {
        for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
		for(par in 1:length(WSIZE)) { # parameter: W_SIZE
			for(MAXGAP.index in 1:length(MAXGAP)){

				print(paste("Parameters READ_LENGTH:",	READLENGTH[[READLENGTH.index]], " W_SIZE:", WSIZE[[par]], " W_PER_C:", WPERC[[par]], " MAX_GAP: ", MAXGAP[[MAXGAP.index]],sep=""))

				runPeakSeq(
					data_directory = data.path,
					work_directory = path.work,
					orig_treatment = treatment.file,
					orig_control = control.file,
					ELAND_PREFIX = treatment.file,
					SGR_PREFIX = treatment.file,
					MAP_FILENAME = map.file,
					INPUT_PREFIX = control.file,
					READ_LENGTH = READLENGTH[[READLENGTH.index]],
					MAX_GAP = MAXGAP[[MAXGAP.index]],
					W_SIZE = WSIZE[[par]],
					W_PER_C = WPERC[[par]],
					OUTPUT = file.path(path.peaks, paste(
						outputName, 
						"_READLENGTH", READLENGTH[[READLENGTH.index]], 
						"_MAXGAP", MAXGAP[[MAXGAP.index]], 
						"_WSIZE", WSIZE[[par]], 
						"_WPERC", WPERC[[par]], 
						"_", sample, sep="")),	
					FINAL = file.path(path.peaks, paste(
						outputName, 
						"_READLENGTH", READLENGTH[[READLENGTH.index]], 
						"_MAXGAP", MAXGAP[[MAXGAP.index]], 
						"_WSIZE", WSIZE[[par]], 
						"_WPERC", WPERC[[par]], 
						"_", sample, "_FINAL", sep="")),	
					LOGFILE = file.path(path.peaks, "log", paste(paste(
						outputName, 
						"_READLENGTH", READLENGTH[[READLENGTH.index]], 
						"_MAXGAP", MAXGAP[[MAXGAP.index]], 
						"_WSIZE", WSIZE[[par]], 
						"_WPERC", WPERC[[par]], 
						"_", sample, "_FINAL", sep="")))	
				)
			}
		} # W_SIZE
	} # READ_LENGTH
}

