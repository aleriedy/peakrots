\name{initialise}
\alias{initialise}
\title{ Initialises peakROTS workflow }
\description{
Initialises peakROTS workflow. Based on the given parameters, this function creates the runtime directory 
structure, workflow description and the settings file. Workflow description contains all jobs with their 
dependencies. The settings file contains all information needed to run the workflow.
}
\arguments{
  \item{detector}{ Peak detection applications, currently supported values are \code{MACS} and \code{PeakSeq}. Default is \code{MACS}.  }
  \item{data.path}{ Path to input data files. Is required. }
  \item{treatment.file}{ Name of the treatment input file. Is required. }
  \item{control.file}{ Name of the control input file. Is required. }
  \item{path.work}{ Path to working directory for storing all runtime information and results. Default is \code{work}. }
  \item{path.bootstrap}{ Path to bootstrap sample directory. Default is \code{<path.work>/bootstrap}. }
  \item{do.bootstrap}{ If true, bootstrap samples are generated, otherwise existing samples are used from \code{path.bootstrap}. Is required. }
  \item{r.command}{ The R executable to use. The R environment must have this library installed. Default is \code{R}. }
  \item{environment.initialiser}{ Shell command that is run before peak detection applications: place platform specific initialisation code here. Default is empty. }
  \item{bootstrap.count}{ Number of individual bootstrap samples. Default is 100. }
  \item{pvalue}{ P-value cutoff for peak detection (MACS). }
  \item{shiftsize}{ List of shiftsizes for peak detection (MACS). }
  \item{tsize}{ List of tsizes for peak detection (MACS) }
  \item{bw}{ List of bw's for peak detection (MACS) }
  \item{nolambda}{ List of values of nolambda parameter for peak detection (MACS) }
  \item{mfold}{ List of mfold values for peak detection (MACS) }
  \item{gsize}{ List of gsize values for peak detection (MACS) }
  \item{READLENGTH}{ List of readlengths for peak detection (PeakSeq) }
  \item{WSIZE}{ List of window sizes for peak detection (PeakSeq) }
  \item{WPERC}{ List of WPERC values for peak detection (PeakSeq) }
  \item{MAXGAP}{ List of MAXGAP values for peak detection (PeakSeq) }
  \item{map.file}{ Mappability file for peak detection (PeakSeq) }
  \item{MAXCHR.defined}{ MAXCHR for peak detection (PeakSeq) }
  \item{preprocess.address}{ PeakSeq preprocessing tool binary path (PeakSeq) }
  \item{peakseq.address}{ PeakSeq binary path (PeakSeq) }
}
\author{ Aleksi Kallio, Teemu D. Laajala }
