#
# Starts running existing peakROTS working directory.
#
run <- function(
	path.work = "work", 
	do.run = do.run.local.silent, 
	polling.timeout.seconds = 1, 
	jobs.running.max = 2,
	verbosity = 0
	) {

	# Load settings
	load.settings(path.work)

	#
	# while there are jobs that are not finished
	#
	while (pendingJobCount() > 0 || runningJobCount() > 0) {	

		if (verbosity > 0) { print("Checking for finished jobs") }

		#
		# check which running jobs have finished
		#
		running.jobs <- readLines(bookkeeping.running)			
		for (job in running.jobs) {
			id <- parse.line(job)["id"]
			lock.file <- file.path(path.work, paste(id, ".finished", sep=""))
			if (file.exists(lock.file)) {
				transition(file.from = bookkeeping.running, file.to = bookkeeping.finished, id = id)
				print(paste("Job ", id, " finished", sep=""))
				unlink(lock.file)
			}
		}

		if (verbosity > 0) { print("Try to start pending jobs") }
		
		#
		# start as many pending jobs as possible and needed
		#
		while (pendingJobCount() > 0 && runningJobCount() < jobs.running.max) {

			# find job that can be run
			if (verbosity > 0) { print("Reading bookkeeping...") }
			pending.jobs <- readLines(bookkeeping.pending)
			finished.job.map <- create.map(readLines(bookkeeping.finished))
			if (verbosity > 0) { print("Bookkeeping read") }

			runnable.job.id <- NULL
			for (job in pending.jobs) {
				id <- parse.line(job)["id"]
				if (verbosity > 1) { print(paste("Checking if job ", id, " can be started", sep="")) }
				dependencies <- parse.line(job)["dependencies"][[1]]
				is.ok <- TRUE
				for (dependency in dependencies) {	
			
					dependency.id <- paste("r", dependency, sep="")
					if (verbosity > 2) { print(paste("Checking if dependency ", dependency.id, " is met", sep="")) }

					# Do lookup from the map
					if (length(finished.job.map) > 0) {
						found <- !is.na(finished.job.map[dependency.id])
					} else {
						found <- FALSE
					}

					if (!found) {
						is.ok <- FALSE
						if (verbosity > 2) { print("Dependency was not met") }
						break
					}
				}
				
				if (is.ok) {
					runnable.job.id <- id
					break
				}
			}
			
			if (!is.null(runnable.job.id)) {
				
				# move pending job to running
				new.job <- transition(file.from = bookkeeping.pending, file.to = bookkeeping.running, runnable.job.id)

				# start the job
				job.command <- paste("echo \"library(peakROTS); dispatch();\" | ", r.command, " --vanilla --args", path.work, new.job["id"], new.job["function"], paste(new.job["args"][[1]], collapse=" "), sep=" ")
				log.file <- file.path(path.log, paste(new.job["id"], ".log", sep=""))
				do.run(job.command, new.job["name"], log.file)
				
				# print status log
				print(paste("Job ", new.job["id"], " started", sep=""))
			
			} else {
				# no of the pending jobs were runnable, so give up and proceed to wait
				break;
			}
		}

		#
		# wait for jobs to finish
		#
		if (pendingJobCount() > 0 || runningJobCount() > 0) {
			if (verbosity > 0) {
				print("Waiting...")
			}

			Sys.sleep(polling.timeout.seconds) # sleep some seconds
		}
	}	
}

jobs.add <- function(file, job.name, script.name, unique.id, dependencies = list(), args = list()) {
	
	# remove lock file for this job, if exists
	lock.file <- file.path(path.work, paste("r", as.character(unique.id), ".finished", sep=""))
	if (file.exists(lock.file)) {
		unlink(lock.file)
	}
	
	# create the line that describes this job
	cat(
			paste("r", as.character(unique.id), sep=""), # unique job id (format: r<id>)
			" ",
			job.name, 	
			" ",
			script.name, 
			" ",
			ifelse(length(dependencies) > 0, paste(paste(dependencies, collapse=" "), " ", sep=""), ""), 			
			"args ",
			ifelse(length(args) > 0, paste(args, collapse=" "), ""), 			
			"\n", 
			file=file, append=TRUE, sep=""
	)
}


create.map <- function(finished.jobs) {
	map <- c()
	for (finished.job in finished.jobs) {
		fields <- strsplit(finished.job, split=" ")[[1]]
		id <- fields[1]
		map[id] <- TRUE;
	}
	return(map)
}


parse.line <- function(line) {
	fields <- strsplit(line, split=" ")[[1]]
	job <- list()
	n.fixed.fields <- 3
	job["id"] <- fields[1]
	job["name"] <- fields[2]
	job["function"] <- fields[3]
	arg.separator <- which(fields == "args")
	job["dependencies"] <- ifelse((n.fixed.fields+1) >= arg.separator, list(), list(fields[(n.fixed.fields+1):(arg.separator-1)]))
	job["args"] <- list(fields[-(1:(arg.separator))])
	return(job)
}

pendingJobCount <- function() {
	return( length(count.fields(bookkeeping.pending)) )
}

runningJobCount <- function() {
	return( length(count.fields(bookkeeping.running)) )
}

transition <- function(id, file.from, file.to) {

	# read from file.from
	lines.from <- readLines(file.from)
		
	# find the line that corresponds to job
	for (i in 1:length(lines.from)) {
		fields <- parse.line(lines.from[i]) 
		if (as.character(id) == as.character(fields["id"])) {
			# found line corresponding to job id 
			line.number <- i
			break;
		}		
	}		
	
	# remove the line from file.from
	job.line <- lines.from[line.number]
	lines.from <- lines.from[-line.number]

	# add the line to file.to and write
	lines.to <- readLines(file.to)	
	lines.to <- c(lines.to, job.line)
	file.create(file.to) # sweep the file
	write.lines(file=file.to, lines.to) # read original lines plus one

	# write file.from
	file.create(file.from) # sweep the file 
	if (length(lines.from) > 0) {
		write.lines(file=file.from, lines.from) # write original lines minus one
	}

	# return parsed version of line that was moved
	return(parse.line(job.line))
}


write.lines <- function(file, lines) {

	con <- file(file)
	writeLines(lines, con = con)
	close(con)
}


do.run.local.verbose <- function(job.command, job.name, log.file) {
	system(job.command, wait=FALSE)	
}

do.run.local <- function(job.command, job.name, log.file) {
	system(paste(job.command, " >> /dev/null", sep=""), wait=FALSE)	
}

do.run.lsf <- function(job.command, job.name, log.file, max.run.time = "24:00") {
	system(paste("bsub -W ", max.run.time, " -J ", job.name, " -u '' -o ", log.file, " -e ", log.file, " '", job.command, "' >> /dev/null", sep=""), wait=FALSE)	
}

do.run.lsf.verbose <- function(job.command, job.name, log.file, max.run.time = "24:00") {
	system(paste("bsub -W ", max.run.time, " -J ", job.name, " -o ", log.file, " -e ", log.file, " '", job.command, "'", sep=""), wait=FALSE)	
}



