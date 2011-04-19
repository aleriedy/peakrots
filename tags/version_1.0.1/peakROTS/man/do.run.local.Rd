\name{do.run.local}
\alias{do.run.local}
\title{ Implements local running of jobs }
\description{
Implements support for running peakROTS jobs locally,
as processes in a single host. Is given as an argument
to function \code{run}. User never calls this function directly.
}
\arguments{
  \item{job.command}{ Job command, passed by peakROTS workflow system }
  \item{job.name}{ Job name, passed by peakROTS workflow system }
  \item{log.file}{ Log file, passed by peakROTS workflow system }
}
\author{ Aleksi Kallio }

