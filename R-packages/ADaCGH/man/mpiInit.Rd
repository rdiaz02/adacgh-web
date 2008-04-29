\name{mpiInit}
\alias{mpiInit}
\title{Initialized a cluster of workstations}
\description{
Use Rmpi and papply, and set a cluster of workstations ready to run ADaCGH.
}
\usage{
mpiInit(wdir = getwd(), minUniverseSize = 15, universeSize = NULL,
        exit_on_fail = FALSE)
}
\arguments{
\item{wdir}{The common ---e.g., NFS mounted--- directory. We need a
  common directory for the graphics, so that they are all found in the
  same location.}
\item{minUniverseSize}{The minimal LAM/MPI universe for the function to
  return successfully. If the function determines that the available
  number of slaves is smaller than minUniverseSize it will fail.}
\item{universeSize}{Alternatively, you can specify the actual size of
  the LAM/MPI universe. If this many slaves cannot be started, the
  funcion will fail.}
\item{exit_on_fail}{If TRUE, kills R session.}
  }
%- maybe also 'usage' for other objects documented here.
\details{
  This function calls Rmpi and related functions to start a cluster of
  workstations, initializes the random number generator, and calls the
  required libraries on the slave nodes.
}
\value{This function is used to create a cluster.}
\references{
Diaz-Uriarte, R. and Rueda, O.M. (2006). ADaCGH: an R package and
web-based application for the analysis of aCGH data. Tech. report.
\url{http://www.ligarto.org/rdiaz/Papers/adacgh.pdf},
\url{http://adacgh.bioinfo.cnio.es}. }
\author{ Ramon Diaz-Uriarte \email{rdiaz02@gmail.com} }
\note{
  We have a somewhat similar function in
  \code{\link[varSelRF]{basicClusterInit.Rd}}, but that is more
  sophisticated.}
\examples{
\dontrun{
mpiInit()
}
}
\keyword{programming}

