
#------------------------------------------------
# FDR-METHODS
#------------------------------------------------

biggerN <- function(x, y) {
 
   ## sort x and y in decreasing order
   x <- sort(x, decreasing = TRUE)		
   y <- sort(y, decreasing = TRUE)		
   ## vector of the positions of (first) matches of the first argument in the second
   a <- match(x, x)				
   ## a logical vector indicating if there is a match or not for its left operand
   b <- x %in% y				
   ## sort c(x,y) in decreasing order
   z <- sort(c(x, y), decreasing = TRUE)	
   ## vector of the positions of (first) matches of the first argument in the second
   d <- match(x, z)				
 
   return(d - a + b)
}

calculateFDR <- function(observed, permuted) {

   ord <- order(observed, decreasing = TRUE)
   a <- observed[ord]
   a.rand <- sort(permuted, decreasing = TRUE)
   n.bigger <- biggerN(a, a.rand)
   A <- vector(length=length(a))
   A[ord] <- rev(cummin(rev(n.bigger/(1:length(a)))))
   ##A[ord] <- n.bigger/(1:length(a))

   return(A)
}



#------------------------------------
# FDR USING THE BOOTSTRAP DATA
#------------------------------------

fdr <- function(selected.shift, selected.tsize, selected.bw, selected.nolambda) {

par1 <- selected.shift  # shift
par2 <- selected.tsize  # tsize
par3 <- selected.bw  # bw
par4 <- selected.nolambda  # nolambda

outputNameO <- "MacspeaksOriginal"
outputNameR <- "MacspeaksRandom"

## Peaks in the original data
peaksO<- parseMACSResults(file.path(path.peaks, paste(outputNameO, 
    "_SHIFT",shiftsize[[par1]], 
    "_TAG",tsize[[par2]], 
    "_BW",bw[[par3]], 
    "_NL",nolambda[[par4]], sep="")))
peaksO[,1]<- sapply(peaksO[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])


## Matrices with rows corresponding to the peaks in the original data and columns corresponding to the randomized detections.
## Here, the results are collected separately for the A- and B-part of the data pairs.
FDR.A<- matrix(nrow=nrow(peaksO), ncol=bootstrap.count)
FDR.B<- matrix(nrow=nrow(peaksO), ncol=bootstrap.count)

for(b in 1:bootstrap.count) {

    ## Peaks in the randomized data for the bootstrap pair b (two randomized datasets for each b).
    peaksAR<- parseMACSResults(file.path(path.peaks, paste(outputNameR, 
        "_SHIFT",shiftsize[[par1]], 
        "_TAG",tsize[[par2]], 
        "_BW",bw[[par3]], 
        "_NL",nolambda[[par4]], 
        "_A", b, sep="")))
    peaksAR[,1]<- sapply(peaksAR[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])

    peaksBR<- parseMACSResults(file.path(path.peaks, paste(outputNameR, 
        "_SHIFT",shiftsize[[par1]], 
        "_TAG",tsize[[par2]], 
        "_BW",bw[[par3]], 
        "_NL",nolambda[[par4]], 
        "_B", b, sep="")))
    peaksBR[,1]<- sapply(peaksBR[,1], function(x) strsplit(x, ".", fixed=TRUE)[[1]][1])

    ## Calculate the FDR for the pair b.
    ## Here, use p-values to order the peaks (note: there are also other possibilities for this).
    FDR.A[,b]<- calculateFDR(observed=peaksO[,"-10xlog10pvalue"], permuted=peaksAR[,"-10xlog10pvalue"])
    FDR.B[,b]<- calculateFDR(observed=peaksO[,"-10xlog10pvalue"], permuted=peaksBR[,"-10xlog10pvalue"])

}

## Calculate FDR as the median over all randomized detections. 
## Set values above 1 to 1.
fdr<- apply(cbind(FDR.A, FDR.B), 1, median, na.rm=TRUE)
fdr[fdr>1]<- 1

save(fdr, file=file.path(path.result, "fdr.Rdata"))


}
