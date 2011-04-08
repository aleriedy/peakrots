bootstrap.dispatch <- function(arglist) {

	print(path.bootstrap)
	# Check environment
	if (!file.exists(path.bootstrap)) {
		dir.create(path.bootstrap)
	}

	sample.number <- as.numeric(arglist[1])

	# TREATMENT FILE
	bootstrap.generate(file.path(data.path, treatment.file), file.path(path.bootstrap, bootstrap.treatment.file), sample.number)

	# CONTROL FILE
	bootstrap.generate(file.path(data.path, control.file), file.path(path.bootstrap, bootstrap.control.file), sample.number)
}

bootstrap.generate <- function(input.file, output.name, sample.number) {

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
		write.table(bsData, file=bsFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	}
}
