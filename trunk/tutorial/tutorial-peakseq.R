
	library(peakROTS)

	initialise(
		detector="PeakSeq",	
		path.work = "work",
		data.path = "data",
		treatment.file = "treatmentdata_example.dat",
		control.file = "controldata_example.dat",
		READLENGTH = list(50),
		MAXGAP = list(100),	
		WSIZE = list(1000000),
		WPERC = list(250),
		map.file = "data/Mapability_CE.txt",
		preprocess.address = "/home/akallio/PeakSeq/Preprocessing/Preprocessing.so",
		peakseq.address = "/home/akallio/PeakSeq/PeakSeq_v1.01/PeakSeq_v1.01.so",
#		bootstrap.count = 10,
		bootstrap.count = 1,
		r.command = "R"
	)

	run(
		path.work = "work", 
		do.run = do.run.local,
		jobs.running.max = 10
	)
	