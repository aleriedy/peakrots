\name{peakROTS-package}
\alias{peakROTS-package}
\alias{peakROTS}
\docType{package}
\title{
ROTS procedure for optimized ChIP-seq peak detection
}
\description{
ROTS procedure for ChIP-seq peak detection allows user to optimize
parameters of common peak detection applications (MACS, PeakSeq). Large optimization
tasks can require significant computational resources. Therefore the package
supports parallel processing either in a single node or multiple nodes (batch processing).
}
\details{
\tabular{ll}{
Package: \tab peakROTS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-10-27\cr
License: \tab GPL (>= 3)\cr
}
The workflow required is computationally heavy, so it cannot be run by a 
single function call. Start by initialising the workflow with \code{initialise}, then 
call \code{run} to start the actual jobs. Depending on the size of the workflow and 
available computational resources running time can be substantial. After
the workflow is finished results can be viewed from the working 
directory (set in \code{path.work} parameter of function \code{initialise}),
from subdirectory \code{results}.
}
\author{
Laura L. Elo, Aleksi Kallio, Teemu D. Laajala

Maintainer: Aleksi Kallio <aleksi.kallio@csc.fi>
}
\references{
Laura L Elo, Aleksi Kallio, Teemu D Laajala, R David Hawkins, Eija Korpelainen & Tero Aittokallio. Optimized detection of transcription factor binding sites in ChIP-seq experiments. Submitted manuscript.

peakROTS web site at: http://www.nic.funet.fi/pub/sci/molbio/peakROTS/
}

