## PeakSeq ROTS, T.D. Laajala, tlaajala@cc.hut.fi
## Summer 2010, Turku University, Math Department


findpeakspeakseq_preprocess.dispatch <- function(arglist) {

	sample <- arglist[5]
	is.original <- (sample == "ORIG") 
	findpeakspeakseq_preprocess(data.path = arglist[1], treatment.file = arglist[2], control.file = arglist[3], outputName = arglist[4], sample = sample, is.original)
}

findpeakspeakseq_preprocess <- function(data.path, treatment.file, control.file, outputName, sample="", is.original) {

	dir.create(file.path(path.peaks, "log"), showWarnings=FALSE, recursive=TRUE)

	    # Creating a temporary folder

	    dir.create(file.path(path.work, "temp"), showWarnings=FALSE, recursive=TRUE)    

	    ## USING GREP TO BREAK DOWN FILES INTO CHROMOSOME SPECIFIC FILES

	    ## 1) CONTROL
	    
	    for(i in 1:22){
		system(paste("more ", file.path(data.path, control.file), " | grep chr", i, ".fa > ", path.work, "/temp/", control.file, "chr", i, ".txt", sep=""))
		print(paste("Using grep on chromosome: ", i, sep=""))
	    }
	    system(paste("more ", file.path(data.path, control.file), " | grep chrX.fa > ", path.work, "/temp/", control.file, "chrX.txt", sep=""))
	    print("Using grep on chromosome: X")
	    system(paste("more ", file.path(data.path, control.file), " | grep chrY.fa > ", path.work, "/temp/", control.file, "chrY.txt", sep=""))
	    print("Using grep on chromosome: Y")
	    system(paste("more ", file.path(data.path, control.file), " | grep chrM.fa > ", path.work, "/temp/", control.file, "chrM.txt", sep=""))
	    print("Using grep on chromosome: M")

	    ## 2) TREATMENT

	    for(i in 1:22){
		system(paste("more ", file.path(data.path, treatment.file), " | grep chr", i, ".fa > ", path.work, "/temp/", treatment.file, "chr", i, ".txt", sep=""))
		print(paste("Using grep on chromosome: ", i, sep=""))
	    }
	    system(paste("more ", file.path(data.path, treatment.file), " | grep chrX.fa > ", path.work, "/temp/", treatment.file, "chrX.txt", sep=""))
	    print("Using grep on chromosome: X")
	    system(paste("more ", file.path(data.path, treatment.file), " | grep chrY.fa > ", path.work, "/temp/", treatment.file, "chrY.txt", sep=""))
	    print("Using grep on chromosome: Y")
	    system(paste("more ", file.path(data.path, treatment.file), " | grep chrM.fa > ", path.work, "/temp/", treatment.file, "chrM.txt", sep=""))
	    print("Using grep on chromosome: M")

	    ## USING PREPROCESSING FUNCTIONS FOR BOTH CONTROL AND TREATMENT
	    # For random peak detection SGR-files have to be generated also for the control

    
	    ## CALLING PREPROCESSING FUNCTIONS

	for(READLENGTH.index in 1:length(READLENGTH)){

	    dir.create(file.path(path.work, "temp", READLENGTH[[READLENGTH.index]]), showWarnings=FALSE, recursive=TRUE)    

	    ## 1) CONTROL

	    print("Preprocessing control...")    
	    print(paste("READ_LENGTH ", READLENGTH[[READLENGTH.index]], sep=""))
	    print(paste("INPUT_FILENAME ", file.path(data.path, control.file), sep=""))
	    print(paste("ELAND_PREFIX ", file.path(path.work, "temp", control.file), sep=""))
	    print(paste("SGR_PREFIX ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file), sep=""))
	    preprocessPeakSeq(READ_LENGTH = READLENGTH[[READLENGTH.index]],
	    		      INPUT_FILENAME= file.path(data.path, control.file),
	    		      ELAND_PREFIX = file.path(path.work, "temp", control.file),
	    		      SGR_PREFIX = file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file))

	    ## 2) TREATMENT

	    print("Preprocessing treatment...")
	    print(paste("READ_LENGTH ", READLENGTH[[READLENGTH.index]], sep=""))
	    print(paste("INPUT_FILENAME ", file.path(data.path, treatment.file), sep=""))
	    print(paste("ELAND_PREFIX ", file.path(path.work, "temp", treatment.file), sep=""))
	    print(paste("SGR_PREFIX ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file), sep=""))

	    preprocessPeakSeq(READ_LENGTH = READLENGTH[[READLENGTH.index]],
			      INPUT_FILENAME= file.path(data.path, treatment.file),
			      ELAND_PREFIX = file.path(path.work, "temp", treatment.file),
			      SGR_PREFIX = file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file))

	}
}

