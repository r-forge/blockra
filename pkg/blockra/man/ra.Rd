% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ra.r
\name{ra}
\alias{ra}
\title{Rearrangement Algorithm}
\usage{
ra(X, f = var, shuffle = T, maxiter = 1000, stalliter = ncol(X),
  abs.tol = 0, rel.tol = 0, f.target = -Inf)
}
\arguments{
\item{f}{function of the rowsums to be minimized}

\item{shuffle}{randomly permute each column of the matrix before rearrangement}

\item{maxiter}{number of maximum iterations}

\item{stalliter}{convergence if no improvement after stalliter iteration}

\item{abs.tol}{abs convergence crit}

\item{rel.tol}{relative convergence criterion}

\item{f.target}{value of the objective function at which absolute converge is reached}
}
\value{
list of rearranged matrix, rowsums (descending), iterations, current objective, iterations of objective, whether converged
}
\description{
Function which performs the traditional RA
}
\author{
Kris Boudt, \email{kris.boudt@vub.ac.be}

Edgars Jakobsons, \email{edgars.jakobsons@math.ethz.ch}

Steven Vanduffel, \email{steven.vanduffel@vub.ac.be}

Kristof Verbeken, \email{kristof.verbeken@vub.ac.be}
}
\references{
\url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2634669}
}
\seealso{
\code{\link{blockra}} for the block rearrangement algorithm

\code{\link{brave}} for the block rearrangement variance equalizer

\code{\link{ra}} for the rearrangement algorithm

\code{\link{blockra}} for the block rearrangement algorithm
}

