
	library(peakROTS)

	initialise(
		path.work = "work",
		data.path = "data",
		treatment.file = "treatmentdata_example.dat",
		control.file = "controldata_example.dat",
#		shiftsize = list(10, 100), 
		shiftsize = list(10),
		tsize = list(25),
		bw = list(300),
		nolambda = list(FALSE),
		pvalue = 10^-2,
#		bootstrap.count = 10,
		bootstrap.count = 1,
		r.command = "R"
	)

	run(
		path.work = "work", 
		do.run = do.run.local,
		jobs.running.max = 10
	)
	