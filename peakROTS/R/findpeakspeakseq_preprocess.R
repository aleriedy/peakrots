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
	    
	    for(i in 1:length(preprocess.chr)){
		#system(paste("more ", file.path(data.path, control.file), " | grep chr", preprocess.chr[i], ".fa > ", file.path(path.work, "temp", control.file), "chr", preprocess.chr[i], ".txt", sep=""))
		system(paste("more ", file.path(data.path, control.file), " | grep ", preprocess.chr[i], ".fa > ", file.path(path.work, "temp", control.file), "chr", preprocess.chr[i], ".txt", sep=""))
		print(paste("Using grep on chromosome: ", preprocess.chr[i], sep=""))
	    }

	    ## 2) TREATMENT

	    for(i in 1:length(preprocess.chr)){
		#system(paste("more ", file.path(data.path, treatment.file), " | grep chr", preprocess.chr[i], ".fa > ", file.path(path.work, "temp", treatment.file), "chr", preprocess.chr[i], ".txt", sep=""))
		system(paste("more ", file.path(data.path, treatment.file), " | grep ", preprocess.chr[i], ".fa > ", file.path(path.work, "temp", treatment.file), "chr", preprocess.chr[i], ".txt", sep=""))
		print(paste("Using grep on chromosome: ", preprocess.chr[i], sep=""))
	    }

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

            # Creating fake SGR files for missing chromosomes so PeakSeq won't crash in peak detection phase due to lack of SGR-file
            # For example C. Elegans can be troublesome without this fake ELAND line.
            for(i in 1:25){
                if(i<23 & length(which(preprocess.chr==i))==0){
                        fakeELAND = paste("LINE\\\tAAAAAAAAAAAAAAAAAAAAAAAAAAA\\\tU0\\\t1\\\t0\\\t0\\\tchr",i,".fa\\\t1000\\\tF\\\t..",sep="")
			fakeSGR = paste("chr",i,"\\\t1000\\\t1",sep="")
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", treatment.file), "chr", i, ".txt", sep=""))
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", control.file), "chr", i, ".txt", sep=""))
			system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file), "chr",i, ".sgr", sep=""))
			system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file), "chr",i, ".sgr", sep=""))
                }
                if(i==23 & length(which(preprocess.chr=="X"))==0){
                        fakeELAND = "LINE\\\tAAAAAAAAAAAAAAAAAAAAAAAAAAA\\\tU0\\\t1\\\t0\\\t0\\\tchrX.fa\\\t1000\\\tF\\\t.."
                        fakeSGR = "chrX\\\t1000\\\t1"
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", treatment.file), "chrX.txt", sep=""))
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", control.file), "chrX.txt", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file), "chrX.sgr", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file), "chrX.sgr", sep="")) 
                }
                if(i==24 & length(which(preprocess.chr=="Y"))==0){
                        fakeELAND = "LINE\\\tAAAAAAAAAAAAAAAAAAAAAAAAAAA\\\tU0\\\t1\\\t0\\\t0\\\tchrY.fa\\\t1000\\\tF\\\t.."
			fakeSGR = "chrY\\\t1000\\\t1"
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", treatment.file), "chrY.txt", sep=""))
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", control.file), "chrY.txt", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file), "chrY.sgr", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file), "chrY.sgr", sep="")) 
                }
                if(i==25 & length(which(preprocess.chr=="M"))==0){
                        fakeELAND = "LINE\\\tAAAAAAAAAAAAAAAAAAAAAAAAAAA\\\tU0\\\t1\\\t0\\\t0\\\tchrM.fa\\\t1000\\\tF\\\t.."
			fakeSGR = "chrM\\\t1000\\\t1"
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", treatment.file), "chrM.txt", sep=""))
                        system(paste("echo ", fakeELAND, " > ", file.path(path.work, "temp", control.file), "chrM.txt", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], treatment.file), "chrM.sgr", sep=""))
                        system(paste("echo ", fakeSGR, " > ", file.path(path.work, "temp", READLENGTH[[READLENGTH.index]], control.file), "chrM.sgr", sep="")) 
                }
            }


	}


}

## The parameters are now parameters of the preprocessing main function and
## this function is called by using the .C-function
preprocessPeakSeq <- function(READ_LENGTH=200, MIN_CNUM=1, MAX_CNUM=23, START=1, STOP=-1,
                              FILENAME_LENGTH=200,
                              INPUT_FILENAME,
                              ELAND_PREFIX,
                              ELAND_SUFFIX=".txt",
                              SGR_PREFIX,
                              SGR_SUFFIX=".sgr"){
    ## Load PeakSeq preprocessing library
    dyn.load(preprocess.address)
                 
  output <- .C("preprocessMain", as.integer(READ_LENGTH), as.integer(MIN_CNUM),
               as.integer(MAX_CNUM), as.integer(START), as.integer(STOP),
               as.integer(FILENAME_LENGTH),
               as.character(INPUT_FILENAME),
               as.character(ELAND_PREFIX),
               as.character(ELAND_SUFFIX),
               as.character(SGR_PREFIX),
               as.character(SGR_SUFFIX))
  return(invisible(output))
}  

