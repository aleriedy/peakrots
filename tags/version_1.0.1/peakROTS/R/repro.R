
repro.dispatch <- function(arglist) {
	# Checking which peak detection method we're using
	print(paste("Running reproducibility for detector: ", detector, sep=""))
	print(paste("Sample number: ", as.numeric(arglist[1]), sep=""))
	print(paste("outputName: ", arglist[2], sep=""))
	print(paste("outputName2: ", arglist[3], sep=""))
	if(detector=="MACS"){
		repro.MACS(sample.number = as.numeric(arglist[1]), outputName = arglist[2], outputName2 = arglist[3])
	}else{
		repro.PeakSeq(sample.number = as.numeric(arglist[1]), outputName = arglist[2], outputName2 = arglist[3])
	}
}

repro.MACS <- function(sample.number, outputName, outputName2) {

	dir.create(file.path(path.repro, "log"), showWarnings=FALSE, recursive=TRUE)

	K<- c(25, c(0.5,1:20)*100, 2500, 3000, 4000, 5000, 7000, c(1:9)*10000)
	sets<- ceiling(c(1:bootstrap.count)/10)

        for(shiftsize.index in 1:length(shiftsize)) { # shiftSize
		for(tsize.index in 1:length(tsize)) { # tSize
	        	for(bw.index in 1:length(bw)) { # bandWidth
	            		for(nolambda.index in 1:length(nolambda)) { # noLambda
	            			for(mfold.index in 1:length(mfold)) { # mfold

								# Matrix for the reproducibility results.			
                        		reproresult<- matrix(nrow=length(K), ncol=5)
                        		colnames(reproresult)<- c("TopK", "p.rand.c500", "p.rand.s500", "e.p.rand.c500", "e.p.rand.s500")
                        		reproresult[,"TopK"]<- K
        
								# Peaks in the first bootstrap dataset.
                        		peaksA<- parseMACSResults(file.path(path.peaks, paste(outputName, 
                            			"_SHIFT",shiftsize[[shiftsize.index]], 
                            			"_TAG",tsize[[tsize.index]], 
                            			"_BW",bw[[bw.index]], 
                            			"_NL",nolambda[[nolambda.index]], 
	        	                        "_MFOLD", mfold[[mfold.index]], 
                            			 paste("_A", sample.number, sep=""), sep="")))
                        		peaksA[,1]<- sapply(peaksA[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])

		                        # Take the central 500bp regions based on peak center.
		                        peaksA.c500<- peaksA
								if (nrow(peaksA) > 0) {
		                     		mm<- floor(apply(cbind(peaksA[,"start"],peaksA[,"end"]),1,median))
		                        	peaksA.c500[,"start"]<- mm - (500-1)
		                        	peaksA.c500[,"end"]<- mm + 500
								}
		                        peaksA.c500<- peaksA.c500[,1:3]

		                        # Take the central 500bp regions based on peak summit.
		                        peaksA.s500<- peaksA
		                        mm<- peaksA[,"start"] + peaksA[,"summit"]
		                        peaksA.s500[,"start"]<- mm - (500-1)
		                        peaksA.s500[,"end"]<- mm + 500
		                        peaksA.s500<- peaksA.s500[,1:3]

								# Random order of the peaks.
		                        set.seed(sample.number) 
		                        randomOrderA<- sample(1:nrow(peaksA), nrow(peaksA))

		                        # Pick 10 bootstrap datasets from the second set of bootstrap datasets.
		                        sel<- which(sets==sets[sample.number])  
		                        for(bb in sel) {
						
		                            	# Peaks in the second bootstrap dataset.
                        			peaksB<- parseMACSResults(file.path(path.peaks, paste(outputName, 
                            				"_SHIFT",shiftsize[[shiftsize.index]], 
	                            			"_TAG",tsize[[tsize.index]], 
	                            			"_BW",bw[[bw.index]], 
	                            			"_NL",nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]], 
	                            			paste("_B", bb, sep=""), sep="")))
        	                		peaksB[,1]<- sapply(peaksB[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])
						
			                        # Take the central 500bp regions based on peak center.
			                        peaksB.c500<- peaksB
									if (nrow(peaksB) > 0) {
										mm<- floor(apply(cbind(peaksB[,"start"],peaksB[,"end"]),1,median))
			                        	peaksB.c500[,"start"]<- mm - (500-1)
			                        	peaksB.c500[,"end"]<- mm + 500
									}
			                        peaksB.c500<- peaksB.c500[,1:3]

			                        # Take the central 500bp regions based on peak summit.
			                        peaksB.s500<- peaksB
			                        mm<- peaksB[,"start"] + peaksB[,"summit"]
			                        peaksB.s500[,"start"]<- mm - (500-1)
			                        peaksB.s500[,"end"]<- mm + 500
			                        peaksB.s500<- peaksB.s500[,1:3]

						 			# Random order of the peaks.
			                        randomOrderB<- sample(1:nrow(peaksB), nrow(peaksB))

                            			# Order the peaks using two different criteria: p-values only or enrichment scores together with p-values.
                            			O.p.rand.A<- order(peaksA[,"-10xlog10pvalue"], randomOrderA, decreasing=TRUE)
                            			O.e.p.rand.A<- order(peaksA[,"fold_enrichment"], peaksA[,"-10xlog10pvalue"], randomOrderA, decreasing=TRUE)
                            			O.p.rand.B<- order(peaksB[,"-10xlog10pvalue"], randomOrderB, decreasing=TRUE)
                            			O.e.p.rand.B<- order(peaksB[,"fold_enrichment"], peaksB[,"-10xlog10pvalue"], randomOrderB, decreasing=TRUE)
                      
                           			# Calculate the reproducibilities at different top list sizes.
                           			for(k in 1:length(K)) {
                                			if(nrow(peaksA)>=K[k] & nrow(peaksB)>=K[k]) {
                                    				reproresult[k,"p.rand.c500"]<- calculateUnionOverlap(as.matrix(peaksA.c500[O.p.rand.A,][1:K[k],]), as.matrix(peaksB.c500[O.p.rand.B,][1:K[k],]))
                                    				reproresult[k,"p.rand.s500"]<- calculateUnionOverlap(as.matrix(peaksA.s500[O.p.rand.A,][1:K[k],]), as.matrix(peaksB.s500[O.p.rand.B,][1:K[k],]))
                                    				reproresult[k,"e.p.rand.c500"]<- calculateUnionOverlap(as.matrix(peaksA.c500[O.e.p.rand.A,][1:K[k],]), as.matrix(peaksB.c500[O.e.p.rand.B,][1:K[k],]))
                                    				reproresult[k,"e.p.rand.s500"]<- calculateUnionOverlap(as.matrix(peaksA.s500[O.e.p.rand.A,][1:K[k],]), as.matrix(peaksB.s500[O.e.p.rand.B,][1:K[k],]))
                                			}
                            			} # k

									# Write result
                	        		write.table(reproresult, file=paste(file.path(path.repro, paste(outputName2, 
                        	    			"_SHIFT",shiftsize[[shiftsize.index]], 
	                        	    		"_TAG", tsize[[tsize.index]], 
	                            			"_BW",bw[[bw.index]], 
        	                    			"_NL",nolambda[[nolambda.index]], 
		        	                        "_MFOLD", mfold[[mfold.index]], 
        	        	            		"_REP", sample.number, 
				                        	"_rep", bb, 
											sep="")), "txt", sep="."), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
		                        } # bb

					} # mfold.index
            			} # noLambda.index
        		} # bw.index
    		} # tSize.index
	} # shiftsize.index
}

repro.PeakSeq <- function(sample.number, outputName, outputName2) {

	dir.create(file.path(path.repro, "log"), showWarnings=FALSE, recursive=TRUE)

	K<- c(25, c(0.5,1:20)*100, 2500, 3000, 4000, 5000, 7000, c(1:9)*10000)
	sets<- ceiling(c(1:bootstrap.count)/10)

        for(READLENGTH.index in 1:length(READLENGTH)) { # parameter: READ_LENGTH
		for(par.index in 1:length(WSIZE)) { # parameter: W_SIZE & W_PER_C (paired)
			for(MAXGAP.index in 1:length(MAXGAP)){

								# Reporting stuff...
					
								print(paste("READLENGTH:", READLENGTH[[READLENGTH.index]], "MAXGAP:", MAXGAP[[MAXGAP.index]],"WSIZE:", WSIZE[[par.index]],"WPERC:",WPERC[[par.index]] ,sep=" ")) 

								# Matrix for the reproducibility results.			
                        		reproresult<- matrix(nrow=length(K), ncol=3)
                        		#colnames(reproresult)<- c("TopK", "p.rand.c500", "p.rand.s500", "e.p.rand.c500", "e.p.rand.s500")
                        		colnames(reproresult)<- c("TopK", "p.rand.c500", "e.p.rand.c500")
                        		reproresult[,"TopK"]<- K
        
								# Peaks in the first bootstrap dataset.
                        		peaksA<- parsePeakSeqResults(file.path(path.peaks, paste(outputName, 
                            			"_READLENGTH",READLENGTH[[READLENGTH.index]], 
                            			"_MAXGAP",MAXGAP[[MAXGAP.index]], 
                            			"_WSIZE",WSIZE[[par.index]], 
                            			"_WPERC",WPERC[[par.index]], 
										paste("_A", sample.number, sep=""),
                            			"_FINAL", sep="")))
                        		peaksA[,1]<- sapply(peaksA[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])

		                        # Take the central 500bp regions based on peak center.
		                        peaksA.c500<- peaksA
								mm<- floor(apply(cbind(peaksA[,"Start"],peaksA[,"Stop"]),1,median))
								peaksA.c500[,"Start"]<- mm - (500-1)
								peaksA.c500[,"Stop"]<- mm + 500
		                        peaksA.c500<- peaksA.c500[,1:3]

		                        # Take the central 500bp regions based on peak summit.
		                        peaksA.s500<- peaksA
		                        # FOR PEAKSEQ THE SUMMIT IS THE SAME AS THE PEAK CENTER
		                        #mm<- peaksA[,"start"] + peaksA[,"summit"]
		                        peaksA.s500[,"Start"]<- mm - (500-1)
		                        peaksA.s500[,"Stop"]<- mm + 500
		                        peaksA.s500<- peaksA.s500[,1:3]

								# Random order of the peaks.
		                        set.seed(sample.number) 
		                        randomOrderA<- sample(1:nrow(peaksA), nrow(peaksA))

		                        # Pick 10 bootstrap datasets from the second set of bootstrap datasets.
		                        sel<- which(sets==sets[sample.number])  
		                        for(bb in sel) {
						
		                            # Peaks in the second bootstrap dataset.
                        			peaksB<- parsePeakSeqResults(file.path(path.peaks, paste(outputName, 
									"_READLENGTH",READLENGTH[[READLENGTH.index]], 
									"_MAXGAP",MAXGAP[[MAXGAP.index]], 
									"_WSIZE",WSIZE[[par.index]], 
									"_WPERC",WPERC[[par.index]], 
									paste("_B", bb, sep=""),
									"_FINAL", sep="")))
        	                		peaksB[,1]<- sapply(peaksB[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])
						
			                        # Take the central 500bp regions based on peak center.
			                        peaksB.c500<- peaksB
									if (nrow(peaksB) > 0) {
			                     		mm<- floor(apply(cbind(peaksB[,"Start"],peaksB[,"Stop"]),1,median))
			                        	peaksB.c500[,"Start"]<- mm - (500-1)
			                        	peaksB.c500[,"Stop"]<- mm + 500
									}
			                        peaksB.c500<- peaksB.c500[,1:3]

			                        # Take the central 500bp regions based on peak summit.
			                        peaksB.s500<- peaksB
			                        # FOR PEAKSEQ THE SUMMIT IS THE SAME AS THE PEAK CENTER
			                        #mm<- peaksB[,"start"] + peaksB[,"summit"]
			                        peaksB.s500[,"Start"]<- mm - (500-1)
			                        peaksB.s500[,"Stop"]<- mm + 500
			                        peaksB.s500<- peaksB.s500[,1:3]

						 			# Random order of the peaks.
			                        randomOrderB<- sample(1:nrow(peaksB), nrow(peaksB))

                            		# Order the peaks using two different criteria: p-values only or enrichment scores together with p-values.
                            		O.p.rand.A<- order(peaksA[,"-10log10q_value"], randomOrderA, decreasing=TRUE)
                            		O.e.p.rand.A<- order(peaksA[,"Enrichment"], peaksA[,"-10log10q_value"], randomOrderA, decreasing=TRUE)
                            		O.p.rand.B<- order(peaksB[,"-10log10q_value"], randomOrderB, decreasing=TRUE)
                            		O.e.p.rand.B<- order(peaksB[,"Enrichment"], peaksB[,"-10log10q_value"], randomOrderB, decreasing=TRUE)
                      
                           			# Calculate the reproducibilities at different top list sizes.
                           			for(k in 1:length(K)) {
                                			if(nrow(peaksA)>=K[k] & nrow(peaksB)>=K[k]) {
                                    				reproresult[k,"p.rand.c500"]<- calculateUnionOverlap(as.matrix(peaksA.c500[O.p.rand.A,][1:K[k],]), as.matrix(peaksB.c500[O.p.rand.B,][1:K[k],]))
                                    				#reproresult[k,"p.rand.s500"]<- calculateUnionOverlap(as.matrix(peaksA.s500[O.p.rand.A,][1:K[k],]), as.matrix(peaksB.s500[O.p.rand.B,][1:K[k],]))
                                    				reproresult[k,"e.p.rand.c500"]<- calculateUnionOverlap(as.matrix(peaksA.c500[O.e.p.rand.A,][1:K[k],]), as.matrix(peaksB.c500[O.e.p.rand.B,][1:K[k],]))
                                    				#reproresult[k,"e.p.rand.s500"]<- calculateUnionOverlap(as.matrix(peaksA.s500[O.e.p.rand.A,][1:K[k],]), as.matrix(peaksB.s500[O.e.p.rand.B,][1:K[k],]))
                                			}
                            			} # k

									# Write result
                	        		write.table(reproresult, file=paste(file.path(path.repro, paste(outputName2, 
										"_READLENGTH",READLENGTH[[READLENGTH.index]], 
										"_MAXGAP",MAXGAP[[MAXGAP.index]], 
										"_WSIZE",WSIZE[[par.index]], 
										"_WPERC",WPERC[[par.index]],
        	        	            	"_REP", sample.number, 
				                    	"_rep", bb, 
										sep="")), "txt", sep="."), quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
		                        } # bb

        		} # MAXGAP.index
    		} # paired W_SIZE & W_PER_C
	} # READLENGTH.index
}

