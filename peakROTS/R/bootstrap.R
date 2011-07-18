bootstrap.dispatch <- function(arglist) {

	print(path.bootstrap)
	# Check environment
	if (!file.exists(path.bootstrap)) {
		dir.create(path.bootstrap)
	}

	sample.number <- as.numeric(arglist[1])

	# CURRENTLY NO PARAMETER IS OFFERED TO GENERATE BOOTSTRAP THE OLD WAY;
	# AS THE MEMORY EFFICIENT VERSION HAS PROVED TO BE ALSO EFFECTIVE IN THE
	# COMPUTATIONAL SPEED SENSE

	# TREATMENT FILE
	bootstrap.generate.memory(file.path(data.path, treatment.file), file.path(path.bootstrap, bootstrap.treatment.file), sample.number)

	# CONTROL FILE
	bootstrap.generate.memory(file.path(data.path, control.file), file.path(path.bootstrap, bootstrap.control.file), sample.number)
}

# File read at whole
bootstrap.generate.speed <- function(input.file, output.name, sample.number) {

	# Read from text file
	# Assume ELAND format
	# 1. Sequence name (derived from file name and line number if format is not Fasta)
	# 2. Sequence
	# 3. Type of match: 
	# 4. Number of exact matches found.
	# 5. Number of 1-error matches found.
	# 6. Number of 2-error matches found.
	# 7. Genome file in which match was found.
	# 8. Position of match (bases in file are numbered starting at 1).
	# 9. Direction of match (F=forward strand, R=reverse).
	# 10. How N characters in read were interpreted: ("."=not applicable, "D"=deletion,
	# 	"I"=insertion).
	# 11. Position and type of first substitution error (e.g. 12A: base 12 was A, not
	# 	whatever is was in read).
	# 12. Position and type of first substitution error, as above. 

	print(paste("Loading ", input.file, "..."))
	data <- read.table(
	input.file, 
	sep="\t", fill=TRUE, header=FALSE, comment.char="",
	col.names=c("Sequence name", "Sequence", "Type of match", "exact matches", "1-error matches", "2-error matches", "Genome file", "Position",
		"Direction", "Interpretation of N characters", "first substitution error", "first substitution error 2"))
	n <- nrow(data)


	## Create bootstrap data files
	for (s in list("A", "B")) {
		bsData <- data[sample(n, n, replace=T), ]
		bsFile <- paste(output.name, ".", s, sample.number, sep="")
		print(paste("Writing file", bsFile, "..."))
		write.table(bsData, file=bsFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
	
	rm(data)
	rm(bsData)
	gc()
}

# Uses less memory (file never read at whole, only single lines)
bootstrap.generate.memory <- function(input.file, output.name, sample.number) {

        line= 0
        print(paste("Counting lines in ", input.file, "..."))
        reading = file(input.file, open="r")

        readLine = readLines(con=reading,n=1)
        while(length(readLine)!=0){
                line=line+1
                readLine = readLines(con=reading,n=1)
        }

        # Total amount of lines known
        n = line
        print(paste("Total amount of lines:",n))
        close(con=reading, type="r")

        ## Create bootstrap data files
        for (s in list("A", "B")) {
                line = 1
                # File and connection
                bsFile <- paste(output.name, ".", s, sample.number, sep="")
                reading = file(input.file, open="r")
                writing = file(bsFile, open="w")

                # Generating bootstrap sample labels
		# Bootstrapping a dataset equal in size to the original dataset
                if(bootstrap.depth==1){
			samples = sample(n, n, replace=T)
		# Bootstrapping a dataset with less tags than original dataset: decreased sequencing depth
		}else if(bootstrap.depth<1){
			samples = sample(n, round(bootstrap.depth*n,0), replace=F)
		# Bootstrapping a dataset with more tags than original dataset: increased sequencing depth
		}else{
			samples = sample(n, round(bootstrap.depth*n,0), replace=T)
		}
                sortedsamples = sort(samples)

                # Writing a line X times (X =0,1,2,3... determined by bootstrap)
                print(paste("Writing file", bsFile, "..."))
                readLine = readLines(con=reading,n=1)
                counter=1
		limit = round(bootstrap.depth*n,0)
                while(counter<limit){
                        if(sortedsamples[counter]==line){
                                writeLines(readLine, con=writing)
                                counter = counter + 1
                        }else{
                        	line = line + 1
                        	readLine = readLines(con=reading,n=1)
                        }
                        if(line%%1000000==0) {print(paste("Line",line, "processed ..."))}
                }
                close(con=writing, type="w")
                close(con=reading, type="r")
        }

}  

