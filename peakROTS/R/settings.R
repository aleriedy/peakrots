settings.file <- "settings.Rdata"

# Load settings to global environment (make them visible everywhere).
load.settings <- function(path.work) {
	load(file.path(path.work, settings.file), envir = globalenv())	
}

# Initialise environment: save setting and create bookkeeping
initialise <- function(
	detector = "MACS",
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

	# Parameters varied to find optimal ROTS

	# Average read length
	#READLENGTH = list(50, 100, 200, 300, 400), # default 200
	# Windowing parameters per chromosome; default assumption that max chromosome size is 250M bp
	#WSIZE = list(250000, 500000, 1000000, 2000000, 5000000), # default 1M
	#WPERC = list(1000, 500, 250, 125, 50), # default 250
	# Maximum gap between peaks
	#MAXGAP = list(100, 200, 300),

	## NEW MORE SOPHISTICATED RUN PARAMETERS TO BETTER OPTIMIZE ROTS
	#READLENGTH = list(50, 100, 150, 200, 250, 300, 350, 400),
	#WSIZE = list(909091,1000000,1111111, 2000000, 5000000),
	#WPERC = list(275,250,225, 125, 50),
	#MAXGAP = list(200),
	
	# Dropping unreliable parameters off
	READLENGTH = list(50, 100, 150, 200, 250, 300, 350, 400),
	WSIZE = list(1000000),
	WPERC = list(250),
	MAXGAP = list (100, 200, 300),	


	# Parameters specific for the run
	
	# Mappability file (different for e.g. human, mice, ...)
	map.file = "/wrk/tlaajala/rots/PeakSeq/Mapability_HG.txt",
	#
	MAXCHR.defined = 23,
	# Preprocess-address
	preprocess.address = "/wrk/tlaajala/rots/run/PeakSeq/Preprocessing/Preprocess.so",
	# PeakSeq-address
	peakseq.address = "/wrk/tlaajala/rots/run/PeakSeq/PeakSeq_v1.01/PeakSeq_v1.01.so"
	
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
		WSIZE,
		WPERC,
		MAXGAP,
		map.file,
		MAXCHR.defined,
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



