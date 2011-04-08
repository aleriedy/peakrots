\name{do.run.lsf}
\alias{do.run.lsf}
\title{ Implements batch execution of jobs (LSF) }
\description{
Implements support for running peakROTS jobs in a LSF
batch processing system. Is given as an argument
to function \code{run}. User never calls this function directly.
}
\arguments{
  \item{job.command}{ Job command, passed by peakROTS workflow system }
  \item{job.name}{ Job name, passed by peakROTS workflow system }
  \item{log.file}{ Log file, passed by peakROTS workflow system }
  \item{max.run.time}{ Maximum running time. Default is \code{24:00}.}
}
\author{ Aleksi Kallio }

