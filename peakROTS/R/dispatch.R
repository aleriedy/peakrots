
dispatch <- function() {

	# parse 3 first command line arguments, because they have special meaning
	args <- commandArgs(trailing=TRUE)

	path.work <- args[1]
	id <- args[2]
	func <- args[3]

	# take rest of the arguments
	args <- args[4:length(args)]

	# load settings
	load.settings(path.work)

	# call the function
	eval(call(paste(func, ".dispatch", sep=""), args))

	# mark this job done
	file.create(file.path(path.work, paste(id, ".finished", sep="")))

}
