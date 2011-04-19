## chipseq testing environment
## version 2009-04-01
## functions for parsing the results files
## ---------------------------------------------------------------------------------------

## If final == FALSE, parse the bootstrap results. Else parse the final results.
## This feature is not yet implemented.


## <name> is the experiment name and it must be same than in the runMACS
## 10xlog10(pvalue) values are sorted decreasingly"
parseMACSResults <- function(name, final=FALSE){
  if( !final ){
    output <- read.table(file=paste(name, "_peaks.xls", sep=""), skip=13, header=TRUE, stringsAsFactors=FALSE)
    ## Choose columns
    output <- output[c("chr","start","end","fold_enrichment", "X.10.log10.pvalue.","tags","summit")]
    ## Fix colnames
    colnames(output)[5] <- "-10xlog10pvalue"
    ## Sort the result according to the -10xlog10(pvalue)
    output <- output[ order(output[,5], decreasing=TRUE), ]
    return(output)
  }
}


## <name> is the location of the output file of the PeakSeq (same as parameter FINAL in
## runPeakSeq)
## q_value values are sorted increasingly
## q_values should be ordered increasingly, but that's not the case if the corresponding
## p-values are equal. It seems the the results are ordered according the p-value.
parsePeakSeqResults <- function(name, final=FALSE){
  if( !final ){
    output <- try(read.table(file=name, header=TRUE, stringsAsFactors=FALSE),silent=TRUE)
    if(class(output)=="try-error"){
    	    # Could not read any lines from input, returning dummy
	    output_new <- matrix(nrow=0, ncol=6)
	    colnames(output_new) = c("Chr", "Start", "Stop","Enrichment", "q_value", "-10log10q_value")
	    return(output_new)
    }else{
	    ## Choose Chr, Start, End, Enrichment and q_value
	    output <- output[c("Chr", "Start", "Stop","Enrichment", "q_value")]
	    # Making a column for -10xlog10q-val
	    output_new <- matrix(nrow=nrow(output), ncol=ncol(output)+1)
	    output_new <- output
	    output_new[,6] <- -10*log10(output[,5])
	    colnames(output_new) = c("Chr", "Start", "Stop","Enrichment", "q_value", "-10log10q_value")
	    return(output_new)
    }
  }
}


## <name> is the analysis path (same as parameter ap in runQuest)
## The result is sorted according to the enrichment intensity of the strongest peak.
parseQuESTResults <- function(name, final=FALSE){
  if( !final ){
    output <- read.table(file=paste(name,"/calls/peak_caller.ChIP.out", sep=""),
                          stringsAsFactors=FALSE)
    ## Find regions
    output <- output[grep("^R-", as.character(output[,1])),]
    ## Parse Start-Stop
    startStop <- as.integer(unlist(strsplit(as.character(output[,3]), "-")))
    startStop <- cbind(startStop[seq(1, length(startStop), by=2)],
                       startStop[seq(2, length(startStop), by=2)])
    output <- cbind(output[,2], startStop, output[,c(5,11)])
    ## Set colnames
    colnames(output) <- c("chr", "Start", "Stop","ChIP","ef")
    return(output)       
  }
}
