results.dispatch <- function(arglist) {

	if(detector=="MACS"){
		results.MACS()
	}else{
		results.PeakSeq()
	}
}

results.MACS <- function() {
	debug <- TRUE
	sets <- ceiling(c(1:bootstrap.count)/10)	

	#--------------------------------------
	# COLLECT THE REPRODUCIBILITY RESULTS
	#--------------------------------------

	## Collect information about the reproducibility results.
	## DD = list of average bootstrap reproducibilities for the different parameter settings
	## DS = list of bootstrap reproducibility standard deviations for the different parameter settings
	## rDD = list of average random reproducibilities for the different parameter settings
	## DZ = list of average bootstrap reproducibility z-scores for the different parameter settings
	## Each element of a list is a matrix corresponding to one parameter setting (one combination of the different parameters) for the peak detection algorithm.
	## The rows of the matrix correspond to different top list sizes, the columns of the matrix correspond to different peak ranking criteria. The first column indicates the top list size.
	#DD<- list()
	#DS<- list()
	#rDD<- list()
	DZ<- list()

	## Names of the parameter settings.
	aa<- list()

	## Index for the parameter setting.
	i<- 1

	outputName <-  "Macsrepro"
	outputName.random <- "MacsreproRandom"

	## For each parameter combination, collect the reproducibility results.
	for(shiftsize.index in 1:length(shiftsize)) { # shiftSize
	    for(tsize.index in 1:length(tsize)) { # tSize
	        for(bw.index in 1:length(bw)) { # bandWidth
	            for(nolambda.index in 1:length(nolambda)) { # noLambda
	            	for(mfold.index in 1:length(mfold)) { # mfold

				if(debug){print(paste("Index for parameter setting: ",i,sep=""))}

        	                ## Reproducibility results for the different bootstrap data pairs and the corresponding randomized data pairs.
                	        apu<- list()
                        	apu.random<- list()
				index <- 1
	                        for(b in 1:bootstrap.count) {
                                    sel <- which(sets==sets[b])
                                    for(bb in sel){
                     
        	                    ## Read the bootstrap reproducibility.
                	            data<- read.table(paste(file.path(path.repro, paste(outputName, 
                        	        "_SHIFT",shiftsize[[shiftsize.index]], 
                                	"_TAG",tsize[[tsize.index]], 
	                                "_BW",bw[[bw.index]],
        	                        "_NL",nolambda[[nolambda.index]], 
        	                        "_MFOLD",mfold[[mfold.index]], 
                	                "_REP",b,
                	                "_rep",bb,
                        	        sep="")), "txt", sep="."), sep="\t", header=TRUE)
	                            apu[[index]]<- data

        	                    ## Read the corresponding random reproducibility.
                	            data.random<- read.table(paste(file.path(path.repro, paste(outputName.random, 
                        	        "_SHIFT",shiftsize[[shiftsize.index]], 
	                                "_TAG",tsize[[tsize.index]], 
        	                        "_BW",bw[[bw.index]], 
                	                "_NL",nolambda[[nolambda.index]], 
        	                        "_MFOLD",mfold[[mfold.index]], 
                        	        "_REP",b,
                	                "_rep",bb,
	                                sep="")), "txt", sep="."), sep="\t", header=TRUE)
        	                    apu.random[[index]]<- data.random

				    index <- index + 1
				  }
                	        }
                            
	                        ## Reproducibility results for the parameter setting indexed i.
        	                ## A.m = average bootstrap reproducibility at different top list sizes
				## A.s = standard deviation of bootstrap reproducibilities at different top list sizes
                        	## rA.m = average random reproducibility at different top list sizes
				## A.z = average bootstrap reproducibility z-score at different top list sizes
				## The rows of each matrix correspond to the different top list sizes, 
				## the columns correspond to different peak ranking criteria used in the reproducibility evaluation;
				## note: currently two unnecessary columns (columns 2 and 4).
				A.m<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
        	                colnames(A.m)<- colnames(apu[[1]])
                	        A.m[,"TopK"]<- apu[[1]][,"TopK"]
                        	A.s<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
	                        colnames(A.s)<- colnames(apu[[1]])
        	                A.s[,"TopK"]<- apu[[1]][,"TopK"]
	                        rA.m<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
        	                colnames(rA.m)<- colnames(apu[[1]])
                	        rA.m[,"TopK"]<- apu[[1]][,"TopK"]
                        	rA.s<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
	                        colnames(rA.s)<- colnames(apu[[1]])
        	                rA.s[,"TopK"]<- apu[[1]][,"TopK"]
                	        A.z<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
                        	colnames(A.z)<- colnames(apu[[1]])
	                        A.z[,"TopK"]<- apu[[1]][,"TopK"]

	                        for(ii in 1:nrow(A.m)) {
        	                    for(jj in 2:ncol(A.m)) {
                	                ## Bootstrap and random reproducibilities for top list size indexed ii and ranking criterion indexed jj.
	                                a<- unlist(sapply(apu, function(x) x[ii,jj]))
        	                        a.random<- unlist(sapply(apu.random, function(x) x[ii,jj]))
                	                ## Number of non-missing values (missing values correspond to cases in which the whole peak list is shorter than the particular top list size).
	                                n1<- sum(!is.na(a))
        	                        n2<- sum(!is.na(a.random))   
                	                ## Require that at least 90% of the values are non-missing.
                        	        if(n1>floor(0.9*bootstrap.count*10) & n2>floor(0.9*bootstrap.count*10)) {
                                	    A.m[ii,jj]<- mean(a, na.rm=TRUE)
	                                    A.s[ii,jj]<- sd(a, na.rm=TRUE)
        	                            rA.m[ii,jj]<- mean(a.random, na.rm=TRUE)
                	                    A.z[ii,jj]<- (A.m[ii,jj]-rA.m[ii,jj])/A.s[ii,jj]
                        	        }
	                           }
        	                }
#	                        DD[[i]]<- A.m
#       	                DS[[i]]<- A.s
#               	        rDD[[i]]<- rA.m
	                        DZ[[i]]<- A.z

	                        aa[[i]]<- paste("_SHIFT",shiftsize[[shiftsize.index]], 
        	                    "_TAG",tsize[[tsize.index]], 
                	            "_BW",bw[[bw.index]], 
                        	    "_NL",nolambda[[nolambda.index]], 
	                            sep="")
        	                i<- i + 1    
			} #nolambda.index    
	            } # par 4
        	} # bw.index
	    } # par 2
	} # shiftsize.index

	#--------------------------------------
	# DETERMINE THE ROTS-RESULTS
	#--------------------------------------

	## Determine the ROTS-result using the ranking criteria indexed jj (currently evaluated jj=3 and jj=5).
	## The ranking criteria are defined in scripts macsreproBootstrap_SHIFT1.R and macsreproBootstrap_random_SHIFT1.R etc.
	#jj<- 5
	jj<- 3

	## Number oftop peaks (kk) considered on the basis of the peak list lengths. Here, require that the peak list with each parameter setting is at least kk. 
	n<- vector(length=length(DZ))
	for(i in 1:length(DZ)) n[i]<- sum(!is.na(unlist(sapply(DZ, function(x) x[i,jj]))))
	kk<- max(DZ[[1]][,"TopK"][n/length(DZ)==1])

	## For each parameter setting, determine the maximum over the different top list sizes and collect the results into a matrix repro.summary.
	sel<- vector(length=length(DZ))
	for(i in 1:length(sel)) {
		m <- which.max(DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj])
		if (length(m) > 0) {
			sel[i] <- m
		} else {
			sel[i] <- -1 # special case: no peaks found
		}
	}
  	repro<- vector(length=length(DZ))
	topks <- vector(length=length(DZ))
	for(i in 1:length(repro)) {
		if (sel[i] != -1) {		
			repro[i] <- DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj][sel[i]]
			topks[i] <- DZ[[1]][,"TopK"][sel[i]]
		} else {
			repro[i] <- 0
			topks[i] <- 0
		}
	}
	repro.summary<- cbind(aa, topks, repro)

	## Choose the most reproducible parameter setting with the two ranking criteria.
	rots<- which.max(repro)
	result <- aa[[rots]]

	dir.create(path.result, showWarnings=FALSE, recursive=TRUE)
	save(result, repro.summary, file=file.path(path.result, "result.Rdata"))
	print(repro.summary)

	## Note: The different lists DD, rDD and DS are collected only for testing different selection criteria for ROTS.
	## They are not needed in the final version.

	return(repro.summary)
}

results.PeakSeq <- function() {
	sets <- ceiling(c(1:bootstrap.count)/10)	

	debug <- FALSE

	#--------------------------------------
	# COLLECT THE REPRODUCIBILITY RESULTS
	#--------------------------------------

	## Collect information about the reproducibility results.
	## DD = list of average bootstrap reproducibilities for the different parameter settings
	## DS = list of bootstrap reproducibility standard deviations for the different parameter settings
	## rDD = list of average random reproducibilities for the different parameter settings
	## DZ = list of average bootstrap reproducibility z-scores for the different parameter settings
	## Each element of a list is a matrix corresponding to one parameter setting (one combination of the different parameters) for the peak detection algorithm.
	## The rows of the matrix correspond to different top list sizes, the columns of the matrix correspond to different peak ranking criteria. The first column indicates the top list size.
	DD<- list()
	#DS<- list()
	#rDD<- list()
	DZ<- list()

	if(debug){print("Setting parameters")}	

	## Names of the parameter settings.
	aa<- list()

	## Index for the parameter setting.
	i<- 1

	outputName <-  "PeakSeqrepro"
	outputName.random <- "PeakSeqreproRandom"

	## For each parameter combination, collect the reproducibility results.
        for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
		for(par.index in 1:length(WSIZE)) { # parameter: W_SIZE & W_PER_C (paired)
			for(MAXGAP.index in 1:length(MAXGAP)){

				if(debug){print("Collecting repro results...")}	


        	                ## Reproducibility results for the different bootstrap data pairs and the corresponding randomized data pairs.
                	        apu<- list()
                        	apu.random<- list()
				index <- 1
	                        for(b in 1:bootstrap.count) {
					sel <- which(sets==sets[b])
					for(bb in sel){
					    # PEAKSEQ: Dropping out columns that would hold summit-information,
					    # as this type of peak data is not available. Only center of peak is used.                     
	
					    # ADDING 39 NA-LINES THE BOTTOM OF REPRO-MATRICES TO GET FULL 75 ROWS 
					    # FOR ALL COMBINATIONS OF PARAMETERS

					    if(debug){print(paste("Reading repro: READLENGTH",READLENGTH[[READLENGTH.index]], "_MAXGAP",MAXGAP[[MAXGAP.index]],"REP",b,"rep",bb,sep=""))}
	        	                    ## Read the bootstrap reproducibility.
        	        	            data<- read.table(paste(file.path(path.repro, paste(outputName, 
						"_READLENGTH",READLENGTH[[READLENGTH.index]], 
						"_MAXGAP",MAXGAP[[MAXGAP.index]], 
						"_WSIZE",WSIZE[[par.index]], 
						"_WPERC",WPERC[[par.index]],
						"_REP",b,
                		                "_rep",bb,
                       		 	        sep="")), "txt", sep="."), sep="\t", header=TRUE)
					    mat2 <- matrix(nrow=39,ncol=3)
					    colnames(mat2)=c("TopK","p.rand.c500","e.p.rand.c500")
				 	    mat2[,1] <- 1000000
	                	            apu[[index]]<-rbind(data,mat2)
					    #if(debug){print(apu[[b]])}

					    if(debug){print(paste("Reading repro random: READLENGTH",READLENGTH[[READLENGTH.index]], "_MAXGAP",MAXGAP[[MAXGAP.index]],"REP",b,"rep",bb,sep=""))}
        		                    ## Read the corresponding random reproducibility.
                		            data.random<- read.table(paste(file.path(path.repro, paste(outputName.random, 
						"_READLENGTH",READLENGTH[[READLENGTH.index]], 
						"_MAXGAP",MAXGAP[[MAXGAP.index]], 
						"_WSIZE",WSIZE[[par.index]], 
						"_WPERC",WPERC[[par.index]],
                   	    	 	        "_REP",b,
                		                "_rep",bb,
	               		                 sep="")), "txt", sep="."), sep="\t", header=TRUE)
					    mat2 <- matrix(nrow=39, ncol=3)
					    colnames(mat2)=c("TopK","p.rand.c500","e.p.rand.c500")
					    mat2[,1] <- 1000000
        	       		            apu.random[[index]]<- rbind(data.random,mat2)
					    #if(debug){print(apu.random[[b]])}

					    index <- index + 1
					}
                	        }

				if(debug){print("Averages, std dev, etc...")}	
                            
	                        ## Reproducibility results for the parameter setting indexed i.
        	                ## A.m = average bootstrap reproducibility at different top list sizes
				## A.s = standard deviation of bootstrap reproducibilities at different top list sizes
                        	## rA.m = average random reproducibility at different top list sizes
				## A.z = average bootstrap reproducibility z-score at different top list sizes
				## The rows of each matrix correspond to the different top list sizes, 
				## the columns correspond to different peak ranking criteria used in the reproducibility evaluation;
				## note: currently two unnecessary columns (columns 2 and 4).
				A.m<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
        	                colnames(A.m)<- colnames(apu[[1]])
                	        A.m[,"TopK"]<- apu[[1]][,"TopK"]
				A.m.s<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
        	                colnames(A.m.s)<- colnames(apu[[1]])
                	        A.m.s[,"TopK"]<- apu[[1]][,"TopK"]
                	        A.s<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
	                        colnames(A.s)<- colnames(apu[[1]])
        	                A.s[,"TopK"]<- apu[[1]][,"TopK"]
	                        rA.m<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
        	                colnames(rA.m)<- colnames(apu[[1]])
                	        rA.m[,"TopK"]<- apu[[1]][,"TopK"]
                        	rA.s<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
	                        colnames(rA.s)<- colnames(apu[[1]])
        	                rA.s[,"TopK"]<- apu[[1]][,"TopK"]
                	        A.z<- matrix(nrow=nrow(apu[[1]]), ncol=ncol(apu[[1]]))
                        	colnames(A.z)<- colnames(apu[[1]])
	                        A.z[,"TopK"]<- apu[[1]][,"TopK"]

				if(debug){print("Loop over whole A.m matrix...")}	

	                        for(ii in 1:nrow(A.m)) {
        	                    for(jj in 2:ncol(A.m)) {
                	                ## Bootstrap and random reproducibilities for top list size indexed ii and ranking criterion indexed jj.
	                                a<- unlist(sapply(apu, function(x) x[ii,jj]))
        	                        a.random<- unlist(sapply(apu.random, function(x) x[ii,jj]))
                	                ## Number of non-missing values (missing values correspond to cases in which the whole peak list is shorter than the particular top list size).
	                                n1<- sum(!is.na(a))
        	                        n2<- sum(!is.na(a.random))   
                	                ## Require that at least 90% of the values are non-missing.
                	            
	
					if(debug){
						print(paste("n1:", n1, "and floor(0.9*bootstrap.count)", floor(0.9*bootstrap.count)))
						print(paste("n2:", n2, "and floor(0.9*bootstrap.count)", floor(0.9*bootstrap.count)))
					}
                	                
					#NOTE! n2 CAN BE PROBLEMATIC! (NRSFOLD, STAT1REP1 and -REP2...)
                        	        if(n1>floor(0.9*bootstrap.count*10) & n2>floor(0.9*bootstrap.count*10)) {
                                	    A.m[ii,jj]<- mean(a, na.rm=TRUE)
	                                    A.s[ii,jj]<- sd(a, na.rm=TRUE)
        	                            rA.m[ii,jj]<- mean(a.random, na.rm=TRUE)
					    # ADDITIONAL CONDITION ADDED BECAUSE OF PROBLEMS WITH RANDOM DATA
					    if(!is.na(rA.m[ii,jj])){
	                	                    A.z[ii,jj]<- (A.m[ii,jj]-rA.m[ii,jj])/A.s[ii,jj]
					    }else{
						    A.z[ii,jj]<- A.m[ii,jj]/A.s[ii,jj]
					    }
					    A.m.s[ii,jj] <- A.m[ii,jj]/A.s[ii,jj]
					   if(debug){
						print(paste("A.m[ii,jj]:", A.m[ii,jj]))
						print(paste("A.s[ii,jj]:", A.s[ii,jj]))
						print(paste("rA.m[ii,jj]:", rA.m[ii,jj]))
						print(paste("A.z[ii,jj]:", A.z[ii,jj]))
					   }
                        	        }
	                           }
        	                }
        	                
				if(debug){print("Results to DZ...")}	
        	                
	                        DD[[i]]<- A.m.s
#       	                DS[[i]]<- A.s
#               	        rDD[[i]]<- rA.m
	                        DZ[[i]]<- A.z

	                        aa[[i]]<- paste(
				    "_READLENGTH",READLENGTH[[READLENGTH.index]], 
			   	    "_MAXGAP",MAXGAP[[MAXGAP.index]], 
				    "_WSIZE",WSIZE[[par.index]], 
				    "_WPERC",WPERC[[par.index]],
	                            sep="")
        	                i<- i + 1    
        		} # MAXGAP.index
    		} # paired W_SIZE & W_PER_C
	} # READLENGTH.index

	#--------------------------------------
	# DETERMINE THE ROTS-RESULTS
	#--------------------------------------

	if(debug){print("Determining results...")}	


	## Determine the ROTS-result using the ranking criteria indexed jj (currently evaluated jj=3 and jj=5).
	## The ranking criteria are defined in scripts macsreproBootstrap_SHIFT1.R and macsreproBootstrap_random_SHIFT1.R etc.
	jj<- 2

	if(debug){print("Number of top peaks (kk) considered on the basis of the peak list lengths")}

	## Number oftop peaks ( kk) considered on the basis of the peak list lengths. Here, require that the peak list with each parameter setting is at least kk. 
	if(debug){print(paste("DZ length:",length(DZ)))}

	# PROBLEMATIC FIELDS!
	# Here length(DZ) is the amount of different combinations (36 for MACS)
	# By chance (or?) the top list sizes TopK happen to also have 36 rows (MACS, produced by repro.R), but
	# for PeakSeq the amount of different parameter combinations aren't 36 but instead 75.
	n<- vector(length=length(DZ))
	n_mean<- vector(length=length(DD))
	for(i in 1:length(DZ)){
		if(debug){
			print(paste("Inside loop, i:", i))
			print(DZ[[i]])
		}
		n[i]<- sum(!is.na(unlist(sapply(DZ, function(x) x[i,jj]))))
		n_mean[i]<- sum(!is.na(unlist(sapply(DD, function(x) x[i,jj]))))
		if(debug){print(paste("n[i]:",n[i]))}
	}
	if(debug){
		print(paste("DZ[[1]]:",DZ[[1]]))
		print(paste("n/length(DZ):",n/length(DZ)))
		print(paste("DZ[[1]][,TopK]:",DZ[[1]][,"TopK"]))
		print(DZ[[1]][,"TopK"][n/length(DZ)==1])
	}
	
	# PROBLEM; NONE WERE FOUND TO BE ==1?? PROBLEMS WITH RANDOM REPRODUCIBILITIES
	kk<- max(DZ[[1]][,"TopK"][n/length(DZ)==1])
	kk_mean<- max(DD[[1]][,"TopK"][n_mean/length(DD)==1])
	
	if(debug){print(paste("kk:",kk))}

	if(debug){print("For each parameter setting determine the maximum over the different top list sizes and collect results")}

	## For each parameter setting, determine the maximum over the different top list sizes and collect the results into a matrix repro.summary.
	sel<- vector(length=length(DZ))
	sel_mean<- vector(length=length(DD))
	for(i in 1:length(sel)){
		if(debug){
			print(paste("Index i:",i))
			print(DZ[[i]][,"TopK"]<=kk)
			print(DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj])
			print(which.max(DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj]))
		}
		if(is.na(DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj][1])){
			sel[i]<- 0
			sel_mean[i] <- 0
		}else{
			sel[i]<- which.max(DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj])
			sel_mean[i] <- which.max(DD[[i]][DD[[i]][,"TopK"]<=kk_mean,jj])
		}
	}
	repro<- vector(length=length(DZ))
	repro_mean<- vector(length=length(DD))
	for(i in 1:length(repro)){
		if(sel[i]==0){
			repro[i]<- 0
		}else{
			repro[i]<- DZ[[i]][DZ[[i]][,"TopK"]<=kk,jj][sel[i]]
			repro_mean[i]<- DD[[i]][DD[[i]][,"TopK"]<=kk_mean,jj][sel_mean[i]]
		}
	}
	repro.summary<- cbind(aa, DZ[[1]][,"TopK"][sel], repro)
	repro_mean.summary<- cbind(aa, DD[[1]][,"TopK"][sel_mean], repro_mean)

	if(debug){print("Choose the most reproducible parameter setting with the two ranking criteria")}

	## Choose the most reproducible parameter setting with the two ranking criteria.
	rots<- which.max(repro)
	rots_mean <- which.max(repro_mean)
	result <- aa[[rots]]
	result_mean <- aa[[rots_mean]]
	

	if(debug){print("Create results")}

	dir.create(path.result, showWarnings=FALSE, recursive=TRUE)
	save(result, repro.summary, result_mean, repro_mean.summary, file=file.path(path.result, "result.Rdata"))

	## Note: The different lists DD, rDD and DS are collected only for testing different selection criteria for ROTS.
	## They are not needed in the final version.

	return(repro.summary)
}
