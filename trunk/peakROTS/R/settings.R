settings.file <- "settings.Rdata"

# Load settings to global environment (make them visible everywhere).
load.settings <- function(path.work) {
	load(file.path(path.work, settings.file), envir = globalenv())	
}

# Initialise environment: save setting and create bookkeeping
initialise <- function(
	detector="MACS",
	data.path,
	treatment.file, 
	control.file, 
	path.work = "work", 
	path.bootstrap = file.path(path.work, "bootstrap"),
	do.bootstrap = TRUE,
	environment.initialiser = "",
	pvalue = 10^-3,
	bootstrap.count = 100,
	r.command = "R",
	
	## MACS ROTS-parameters
	shiftsize = list(NULL, 1, 50, 100, 150, 200), 	# Shift size (shiftsize): 10, 100, 200, 500 (--nomodel=True), shift size must be > 0
	tsize = list(25),	 	# Tag size (tsize): 25, 50, 100, 200, 500 (default 25)
	bw = list(100, 300, 500), 	# Bandwidth (bw): (default 300)
	nolambda = list(FALSE, TRUE), 	# Use of local lambda (nolambda): With this flag on, MACS will use the background lambda as local lambda
	mfold = list(32),
	gsize = 2700000000,       # Effective genome size, default:2.7e+9
	
	## PeakSeq ROTS-parameters

	# Assuming default WSIZE & WPERC
	READLENGTH = list(50, 100, 150, 200, 250, 300, 350, 400),
	MAXGAP = list(100, 200, 300),	
	WSIZE = list(1000000),
	WPERC = list(250),


	# Parameters specific for the run
	# User has to set these to suit the analysis to be run
	
	# Mappability file (different for e.g. human, mice, ...)
	# For example:
	# map.file = "/wrk/tlaajala/rots/PeakSeq/Mapability_HG.txt",
	map.file,
	
	# Max Chromosome to be included in study
	MAXCHR.defined = 23,
	
	# Chromosomes to go through in preprocessing: 
	# for example C. Elegans coding is I, II, III, IV, V, X, M instead
	# of the style chromosome 1,2,..., X, Y, M for human
	preprocess.chr = as.character(c(1:22, "X", "Y", "M")),
	# For example C. Elegans:
	# preprocess.chr = as.character(c("I", "II", "III", "IV", "V", "X", "M"))
	# or 
	# preprocess.chr = as.character(c(1:5, "X", "M"))
	# to match the columns with "chrIV.fa" or "chr4.fa" in the ELAND file.
	

	# Preprocess-address:
	# The address to the shared library "Preprocess.so":
	# e.g. preprocess.address = "/wrk/tlaajala/rots/run/PeakSeq/Preprocess/Preprocess.so"
	preprocess.address,
	
	# PeakSeq-address:
	# The address to the shared library "PeakSeq_v1.01.so":
	# e.g. peakseq.address = "/wrk/tlaajala/rots/run/PeakSeq/PeakSeq_v1.01/PeakSeq_v1.01.so"
	peakseq.address
	
	) {

	# Directory structure and file names	
	path.peaks <- file.path(path.work, "peaks")
	path.repro <- file.path(path.work, "repro")
	path.result <- file.path(path.work, "result")
	path.log <- file.path(path.work, "log")
		
	bookkeeping.pending = file.path(path.work, "pending-jobs.txt")
	bookkeeping.running = file.path(path.work, "running-jobs.txt")
	bookkeeping.finished = file.path(path.work, "finished-jobs.txt")
	
	bootstrap.treatment.file <- "treatmentdata.dat"
	bootstrap.control.file <- "controldata.dat"

	# File formats
	format <- "ELAND"
	
	# Create shared parts of the directory structure
	for (d in c(path.work, path.log)) {
		if (!file.exists(d)) {
			dir.create(d)
		}
	}

	# Save settings to file
	save(file=file.path(path.work, settings.file), 
		detector,
		path.work, 
		path.bootstrap, 
		data.path,
		treatment.file, 
		control.file, 
		path.peaks,
		path.repro,
		path.result,
		path.log,
		format,
		bootstrap.count,
		r.command,
		pvalue,
		shiftsize,
		tsize,
		bw,
		nolambda,
		mfold,
		READLENGTH,
		MAXGAP,
		WSIZE,
		WPERC,
		map.file,
		MAXCHR.defined,
		preprocess.chr,
		preprocess.address,
		peakseq.address,
		bookkeeping.pending,
		bookkeeping.running,
		bookkeeping.finished,
		bootstrap.treatment.file,
		bootstrap.control.file,
		environment.initialiser,
		gsize
)

	# Make settings available in the caller's environment
	load.settings(path.work)

	# Fill pending jobs
	file.create(bookkeeping.pending)
	rots.create.jobs(bookkeeping.pending, do.bootstrap = do.bootstrap)

	# Initialise rest of bookkeeping
	file.create(bookkeeping.running)
	file.create(bookkeeping.finished)
}



