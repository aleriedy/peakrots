## PeakSeq ROTS, T.D. Laajala, tlaajala@cc.hut.fi
## Summer 2010, Turku University, Math Department


# 3 Stages of PeakSeq run:
#
# 1: Use grep to break down original .dat-file into chromosome-specific files
# 2: Use Preprocess-function on the chromosome-specific files to create SGR files
# 3: Finally call the PeakSeqMain function to do peak detection 


runPeakSeq <-
  function(
  	   data_directory,
  	   work_directory,
  	   orig_control,
  	   orig_treatment,
           ## The first chromosome used. (23 is X, 24 is Y, 25 is M.)
           MIN_CNUM=1,
           ## The last chromosome used. (23 is X, 24 is Y, 25 is M.)
           MAX_CNUM=MAXCHR.defined,
           ## Number of simulations per window to estimate FDR.
           N_SIMS=10,
           ## Average length of DNA fragments.
           READ_LENGTH=200,
           ## The maximum gap allowed between peaks for them to be merged together.
           ## Hits that have greater separations aren't merged.
           MAX_GAP=200,
           ## Required false discovery rate.
           MIN_FDR=0.05,
           ## Maximum number of reads using the same starting nulceotide.
           MAX_COUNT=100,
           ## Window size for scoring (in nucleotides).
           W_SIZE=1000000,
           ## Number of windows per chromosome
           W_PER_C=250,
           ## Bin size for doing linear regression.
           BIN_SIZE=10000,
           ## Bin size for doing linear regression for mitochondria.
           BIN_SIZE_M=1000,
           ## The amount on each side that regions are extended when extended regions
           ## are used.
           EXTENDED_REGION_SIZE=2000,
           ## The threshold pvalue for a peak to be outputted to the final file.
           PVAL_THRESH=0.2,
           ## The largest region for which the extended region is evaluated. For peaks
           ## larger than this, only the peak is evaluated for the control input file.
           MAX_REG_EXT=1000,
           ## The following are used to construct filenames. A string of the format chr1
           ## to chrM is placed between the PREFIX and SUFFIX of each filename.
           ## The sample Eland files.
           ELAND_PREFIX,
           ELAND_SUFFIX=".txt",
           ## The sample sgr files.
           SGR_PREFIX,
           SGR_SUFFIX=".sgr",
           ## The file with the information about the mappability fractions.
           MAP_FILENAME,
           ## The control input Eland files.
           INPUT_PREFIX,
           INPUT_SUFFIX=".txt",
           ## The location of the intermediate output file with peaks sorted by
           ## position.
           OUTPUT,
           ## The location of the final output file with peaks sorted by q-value
           ## (p-value after BH correction for multiple hypothesis testing).
           FINAL,
           ## Parameters and calculations based on parameters that are not to be
           ## changed. Minimum allowed threshold.
           MIN_THRESH=2,
           ## Maximum allowed threshold.
           MAX_THRESH=100,
           ## The total number of thresholds.
           N_THRESHES=(MAX_THRESH - MIN_THRESH + 1),
           ## The total number of chromosomes used in the program.
           N_CHRS=(MAX_CNUM - MIN_CNUM + 1),
           ## The number of bins on a single chromosome.
           N_BINS=W_SIZE * W_PER_C / BIN_SIZE,
           ## The maximum number of characters in a filename.
           FILENAME_LENGTH=200,
           ## The score of a start position of a read used in determining the aggregate
           ## eight of a position.
           START=1,
           ## The score of a stop position of a read used in determining the aggregate
           ## height of a position.
           STOP=-1,
           LOGFILE="/dev/null"){

    ## THE BELOW PATHS HAVE TO BE PREDEFINED SO THAT PEAKSEQ SHARED LIBRARIES ARE AVAILABLE
    ## Load PeakSeq preprocessing library
    dyn.load(preprocess.address)
    ## Load PeakSeq
    dyn.load(peakseq.address)
    
    ## CALLING ACTUAL PEAK DETECTION FUNCTION
    print("Commencing actual peak detection...")    
    output <- .C("PeakSeqMain",
                 as.integer(MIN_CNUM), as.integer(MAX_CNUM), as.integer(N_SIMS),
                 as.integer(READ_LENGTH), as.integer(MAX_GAP), as.double(MIN_FDR),
                 as.integer(MAX_COUNT), as.integer(W_SIZE), as.integer(W_PER_C),
                 as.integer(BIN_SIZE), as.integer(BIN_SIZE_M),
                 as.integer(EXTENDED_REGION_SIZE), as.double(PVAL_THRESH),
                 as.integer(MAX_REG_EXT), as.character(file.path(work_directory,"temp",ELAND_PREFIX)),
                 as.character(ELAND_SUFFIX), as.character(file.path(work_directory,"temp",READ_LENGTH,SGR_PREFIX)),
                 as.character(SGR_SUFFIX), as.character(MAP_FILENAME),
                 as.character(file.path(work_directory,"temp",INPUT_PREFIX)), as.character(INPUT_SUFFIX),
                 as.character(OUTPUT), as.character(FINAL),
                 as.integer(MIN_THRESH), as.integer(MAX_THRESH),
                 as.integer(FILENAME_LENGTH), as.integer(START), as.integer(STOP),
                 as.character(LOGFILE))

    ## TODO!! CLEAN UP - REMOVING TEMPORAL FILES, folder /temp/ in working folder
    #print("Clean up, removing temporal files...")    


    return(invisible(output))
  }


## The parameters are now parameters of the preprocessing main function and
## this function is called by using the .C-function
preprocessPeakSeq <- function(READ_LENGTH=200, MIN_CNUM=1, MAX_CNUM=23, START=1, STOP=-1,
                              FILENAME_LENGTH=200,
                              INPUT_FILENAME,
                              ELAND_PREFIX,
                              ELAND_SUFFIX=".txt",
                              SGR_PREFIX,
                              SGR_SUFFIX=".sgr"){
  
  output <- .C("preprocessMain", as.integer(READ_LENGTH), as.integer(MIN_CNUM),
               as.integer(MAX_CNUM), as.integer(START), as.integer(STOP),
               as.integer(FILENAME_LENGTH),
               as.character(INPUT_FILENAME),
               as.character(ELAND_PREFIX),
               as.character(ELAND_SUFFIX),
               as.character(SGR_PREFIX),
               as.character(SGR_SUFFIX))
  return(invisible(output))
}


testPeakSeq <- function(){
	runPeakSeq(data_directory = "/wrk/lelo/chipseq-rots-data/data/",
		   work_directory = "/wrk/tlaajala/rots/",
		   orig_control = "controldata_FoxA1.dat",
		   orig_treatment = "treatmentdata_FoxA1.dat",
		   ELAND_PREFIX = "treatmentdata_FoxA1.",
		   SGR_PREFIX = "treatmentdata_FoxA1.",
		   MAP_FILENAME = "/wrk/tlaajala/rots/PeakSeq/Mapability_HG.txt",
		   INPUT_PREFIX = "controldata_FoxA1.",
		   OUTPUT="/wrk/tlaajala/rots/FoxA1_test.txt",
		   FINAL="/wrk/tlaajala/rots/FoxA1_test_final.txt")

}
