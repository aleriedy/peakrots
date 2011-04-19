convert.to.eland <- function(input, output) {


	## Read from text file and save
    	## ELAND format
    	## 1. Sequence name
    	## 2. Sequence
    	## 3. Type of match: 
    	## 4. Number of exact matches found.
    	## 5. Number of 1-error matches found.
    	## 6. Number of 2-error matches found.
    	## 7. Genome file in which match was found.
    	## 8. Position of match (bases in file are numbered starting at 1).
    	## 9. Direction of match (F=forward strand, R=reverse).
    	## 10. How N characters in read were interpreted: ("."=not applicable, "D"=deletion,
    	##     "I"=insertion).
    	## 11. Position and type of first substitution error (e.g. 12A: base 12 was A, not
    	##     whatever is was in read).
    	## 12. Position and type of first substitution error, as above. 

	data <- read.table(input, sep="\t", fill=TRUE, header=FALSE, comment.char="", stringsAsFactors=FALSE)

	a<- data[,3]
	a[a=="+"]<- "F"
	a[a=="-"]<- "R"

	data<- cbind(
	    paste("LINE", 1:nrow(data), sep=""), 
	    rep("AAAAAAAAAAAAAAAAAAAAAAAAA", nrow(data)),
	    rep("U0", nrow(data)),
	    rep(1, nrow(data)),
	    rep(0, nrow(data)),
	    rep(0, nrow(data)),
	    paste(data[,1], ".fa", sep=""),
	    data[,2],
	    a,
	    rep("..", nrow(data)))
	write.table(data, file=output, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

peakamounts <- function(path.work="/wrk/tlaajala/rots/new/work-foxa1-peakseq-new/", run=c(TRUE,FALSE,FALSE)){

	# AMOUNT OF DIFFERENT PARAMETER SETTINGS FOR PEAKSEQ:
	ORIG_amount <- 24
	RANDOM_amount <- 15000
	BOOTST_amount <- 15000
	
	b <- 1:100
	s <- c("A", "B")
	index <- 1

	cols <- 0
	if(run[1]){cols <- cols + ORIG_amount}
	if(run[2]){cols <- cols + BOOTST_amount}
	if(run[3]){cols <- cols + RANDOM_amount}

	peak_amounts <- matrix(nrow=cols, ncol=2)
	#rownames(peak_amounts) <- c(rep("Original", 75), rep("Bootstrap",15000), rep("Random",15000))

	if(run[1]==TRUE){
		# ORIG peaks
		prefix <- "PeakSeqpeaksOriginal"
		for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
			for(par.index in 1:length(WSIZE)) { # parameter: W_SIZE & W_PER_C (paired)
				for(MAXGAP.index in 1:length(MAXGAP)){
								filename <- paste(
								prefix,
								"_READLENGTH",READLENGTH[[READLENGTH.index]], 
								"_MAXGAP",MAXGAP[[MAXGAP.index]], 
								"_WSIZE",WSIZE[[par.index]], 
								"_WPERC",WPERC[[par.index]], 
								"_ORIG",
								"_FINAL", sep="")		
								temp <- read.table(file=paste(path.work, "peaks/", filename, sep=""), header=TRUE)
								peak_amounts[index,1] <- length(temp[,1])
								peak_amounts[index,2] <- filename
								#print(peak_amounts[index,])
								index <- index +1
				}
			}
		}
	}
	if(run[2]==TRUE){
		# Bootstrap peaks
		prefix <- "PeakSeqpeaks"
		for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
			for(par.index in 1:length(WSIZE)) { # parameter: W_SIZE & W_PER_C (paired)
				for(MAXGAP.index in 1:length(MAXGAP)){
					for(bb in 1:length(b)){
						for(ss in 1:length(s)){
								filename <- paste(
								prefix,
								"_READLENGTH",READLENGTH[[READLENGTH.index]], 
								"_MAXGAP",MAXGAP[[MAXGAP.index]], 
								"_WSIZE",WSIZE[[par.index]], 
								"_WPERC",WPERC[[par.index]], 
								paste("_", s[ss], b[bb], sep=""),
								"_FINAL", sep="")
								temp <- read.table(file=paste(path.work, "peaks/", filename, sep=""), header=TRUE)
								peak_amounts[index,1] <- length(temp[,1])
								peak_amounts[index,2] <- filename
								#print(peak_amounts[index,])
								index <- index +1                            			
						}
					}	
				}
			}
		}
	}
	if(run[3]==TRUE){
		# Random peaks
		prefix <- "PeakSeqpeaksRandom"
		for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
			for(par.index in 1:length(WSIZE)) { # parameter: W_SIZE & W_PER_C (paired)
				for(MAXGAP.index in 1:length(MAXGAP)){
					for(bb in 1:length(b)){
						for(ss in 1:length(s)){
								filename <- paste(
								prefix,
								"_READLENGTH",READLENGTH[[READLENGTH.index]], 
								"_MAXGAP",MAXGAP[[MAXGAP.index]], 
								"_WSIZE",WSIZE[[par.index]], 
								"_WPERC",WPERC[[par.index]], 
								paste("_", s[ss], b[bb], sep=""),
								"_FINAL", sep="")
								temp <- read.table(file=paste(path.work, "peaks/", filename, sep=""), header=TRUE)
								peak_amounts[index,1] <- length(temp[,1])
								peak_amounts[index,2] <- filename
								#print(peak_amounts[index,])
								index <- index +1                            			
						}
					}	
				}
			}
		}
		#print(paste("Min amount:", min(peak_amounts[,1])))
		#print(paste("Max amount:", max(peak_amounts[,1])))
		#print(paste("Avg amount:", mean(peak_amounts[,1])))
		#print(paste("Med amount:", median(peak_amounts[,1])))
	}
	write.table(peak_amounts, file=paste(path.work, "/peak_amounts.txt",sep=""), quote=FALSE, sep="\t")
}


