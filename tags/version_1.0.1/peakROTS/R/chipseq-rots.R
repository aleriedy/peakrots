
#
# The ChIP-Seq ROTS high-level execution plan
# 

rots.create.jobs <- function(file, do.bootstrap) {
	
	unique.id <- 1
	
	# Bootstrap sample generation
	if (do.bootstrap) {
		for (b in 1:bootstrap.count) {
			jobs.add(
				file = file, 
				job.name = paste("bootstrap", as.character(b), sep=""), 
				script.name = "bootstrap", 
				unique.id = unique.id, 
				args = as.character(b)
			)
			unique.id <- unique.id + 1
		}
	}

	# Peak finding for the original data
	if(detector=="MACS"){
		jobs.add(file = file, job.name = "findpeaks_orig", script.name = "findpeaksmacs", unique.id = unique.id, args = list(file.path(data.path, treatment.file), file.path(data.path, control.file), "MacspeaksOriginal", "ORIG"))
	} else{
		# Call preprocessing for the very first file, after this file the preprocessing tasks are to be queued
		
		jobs.add(file = file, job.name = "peakseq_preprocess", script.name = "findpeakspeakseq_preprocess", unique.id = unique.id, args = list(data.path, treatment.file, control.file, "PeakSeqpeaksOriginal", "ORIG"))
		unique.id.first <- unique.id
		unique.id <- unique.id + 1
		jobs.add(file = file, job.name = "findpeaks_orig", script.name = "findpeakspeakseq", unique.id = unique.id, dependencies=unique.id.first, args = list(data.path, treatment.file, control.file, "PeakSeqpeaksOriginal", "ORIG"))
	}
	unique.id.findpeaks.original <- unique.id
	unique.id <- unique.id + 1


	# If we're using PeakSeq adding preprocessing tasks for Bootstrap/Random to queue
	id.before.findpeaks <- unique.id
	if(detector=="PeakSeq"){
		for (b in 1:bootstrap.count) {
			for (s in list("A", "B")) {
				args = list(			
					path.bootstrap,
					paste(bootstrap.treatment.file, ".", s, b, sep=""), 
					paste(bootstrap.control.file, ".", s, b, sep=""), 
					"PeakSeqpeaks", 
					paste(s, b, sep="")			
				)

				#jobs.add(file = file, job.name = "peakseq_preprocess", script.name = "findpeakspeakseq_preprocess", unique.id = unique.id, args)

				jobs.add(
					file = file, 
					job.name = paste("peakseq_preprocess", s, b, sep=""), 
					script.name = "findpeakspeakseq_preprocess",
					unique.id = unique.id, 
					# Initial peak detection should be done first
					dependencies = 1:2,
					args = args
				)

				unique.id <- unique.id + 1			
			}
		}
	}

	id.after.preprocessing <- unique.id

	## DEFINING DEPENDENCIES IN A NEW WAY
	new.dependencies = 1:(id.after.preprocessing-1)
	
	# Peak finding for bootstrap and randomized data
	for (b in 1:bootstrap.count) {
		for (s in list("A", "B")) {

			# Bootstrap
			if(detector=="MACS"){
				script.name = "findpeaksmacs"
				args = list(
					paste(file.path(path.bootstrap, bootstrap.treatment.file), ".", s, b, sep=""), 
					paste(file.path(path.bootstrap, bootstrap.control.file), ".", s, b, sep=""), 
					"Macspeaks", 
					paste(s, b, sep="")
				)
			}else{
				script.name = "findpeakspeakseq"
				args = list(
					path.bootstrap,
					paste(bootstrap.treatment.file, ".", s, b, sep=""), 
					paste(bootstrap.control.file, ".", s, b, sep=""), 
					"PeakSeqpeaks", 
					paste(s, b, sep="")
				)
			}

			# Bootstrap
			jobs.add(
				file = file, 
				job.name = paste("findpeaks_b", s, b, sep=""), 
				script.name = script.name,
				unique.id = unique.id, 
				dependencies = new.dependencies,
				args = args
			)
			unique.id <- unique.id + 1

			# Randomized
			if(detector=="MACS"){
				script.name = "findpeaksmacs"
				args = list(
					paste(file.path(path.bootstrap, "controldata.dat"), ".", s, b, sep=""), 
					paste(file.path(path.bootstrap, "treatmentdata.dat"), ".", s, b, sep=""), 
					"MacspeaksRandom", 
					paste(s, b, sep="")
				)
			} else{
				script.name = "findpeakspeakseq"
				args = list(
					path.bootstrap,
					paste("controldata.dat", ".", s, b, sep=""), 
					paste("treatmentdata.dat", ".", s, b, sep=""), 
					"PeakSeqpeaksRandom", 
					paste(s, b, sep="")
				)
			}
			# Randomized (call MACS/PeakSeq with treatment/control swapped)
			jobs.add(
				file = file, 
				job.name = paste("findpeaks_r", s, b, sep=""), 
				script.name = script.name,
				unique.id = unique.id,
				dependencies = new.dependencies,
				args = args
			)
			unique.id <- unique.id + 1
		}
	}
	findpeaks.ids <- id.before.findpeaks:(unique.id-1)

	# Peak reproducibility for bootstrap and randomized data
	id.before.repro <- unique.id
	for (b in 1:bootstrap.count) {
		# Bootstrap
		if(detector=="MACS"){
			args = list(
				as.character(b), 
				"Macspeaks",
				"Macsrepro"
			)
		}else{
			args = list(
				as.character(b), 
				"PeakSeqpeaks",
				"PeakSeqrepro"
			)
		}
		jobs.add(
			file = file, 
			job.name = paste("repro_b",  b, sep=""), 
			script.name = "repro", 
			unique.id = unique.id, 
			dependencies = findpeaks.ids, # repro is free to use any of the peak results
			args
		)
		unique.id <- unique.id + 1

		# Randomized
		if(detector=="MACS"){
			args = list(
				as.character(b), 
				"MacspeaksRandom",
				"MacsreproRandom"
			)
		}else{
			args = list(
				as.character(b), 
				"PeakSeqpeaksRandom",
				"PeakSeqreproRandom"
			)
		}
		jobs.add(
			file = file, 
			job.name = paste("repro_r", b, sep=""), 
			script.name = "repro", 
			unique.id = unique.id, 
			dependencies = findpeaks.ids, # repro is free to use any of the peak results
			args
		)
		unique.id <- unique.id + 1
	}
	repro.ids <- id.before.repro:(unique.id-1)

	# Collect results
	jobs.add(
		file = file, 
		job.name = "results",
		script.name = "results", 
		unique.id = unique.id, 
		dependencies = repro.ids,
		args = list()
	)
	unique.id <- unique.id + 1
}




