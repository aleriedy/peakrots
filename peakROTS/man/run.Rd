\name{run}
\alias{run}
\title{ Runs a previously initialised peakROTS workflow}
\description{
Runs (starts or resumes) a peakROTS workflow that has been initialised by
a call to function \code{initialise}. Actual implementation 
depends on the \code{do.run} function that is passed.
This function can be interrupted and it will 
continue when called for the next time.
}
\arguments{
  \item{path.work}{ Path to peakROTS workflow, i.e., working directory  }
  \item{do.run}{ The function that is called for each job to actually run it }
  \item{polling.timeout.seconds}{ How often this function polls running jobs for state changes }
  \item{jobs.running.max}{ Maximum number of simultaneous jobs }
  \item{verbosity}{ Level of verbosity for this function. Possible values are 0, 1, 2 and 3.  }
}
\details{
State of the workflow (i.e., states of the individual jobs) is kept on disk. This function
can be interrupted and called again later, as it resumes the state from
the disk. However there should be only one instance of this function
running against the same workflow at the same time. \code{do.run} function can be 
changed between calls to this function.
}
\author{ Aleksi Kallio }

