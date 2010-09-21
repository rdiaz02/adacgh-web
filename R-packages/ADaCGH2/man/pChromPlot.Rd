\name{pChromPlot}
\alias{pChromPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Segment plots for aCGH}


\description{ Segment plots (by chromosome) for aCGH segmentation
  results. Internal calls are parallelized for increased speed and we
  use ff objets to allow the handling of very large objects. The output
  produces files (and associated png) with an image m} }


\usage{
pChromPlot(outRDataName, cghRDataName, chromRDataName,
           probenamesRDataName,
           posRDataName = NULL,
           imgheight = 500,
           pixels.point = 3,
           pch = 20,
           colors = c("orange", "red", "green", "blue", "black"),
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outRDataName}{
    The Rdata file name that contains the
    \code{\link[ff]{ffdf}} with the results from the segmentation as
    carried out by any of the \code{\link{pSegment}} functions.
  }


  \item{cghRDataName}{The Rdata file name that contains the
    \code{\link[ff]{ffdf}} with the aCGH data. This file can be created
    using \code{\link[ff]{as.ffdf}} with a data frame with genes
    (probes) in rows and subjects or arrays in columns. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files.}
  
  
  \item{chromRDataName}{The RData file name with the ff (short integer)
    vector with the chromosome indicator. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files. }
  
 

  \item{probenamesRDataName}{The RData file name with the 
    vector with the probe names. Function
    \code{\link{inputDataToADaCGHData}} produces these type of
    files. Note that this is not an ff file.

}

 \item{posRDataName}{The RData file name with the ff double vector with
    the location (e.g., position in kbases) of each probe in
    the chromosome. Function
    \code{\link{inputDataToADaCGHData}} produces these type of files.
    Used for the spacing in the plots. If NULL, the x-axis goes from
    1:number of probes in that chromosome.
  }


  \item{imgheight}{
    Height of png image. See \code{\link{png}}.
  }
  
  \item{pixels.point}{
    Approximate number of pixels that each point takes; this determines
    also final figure size. With many probes per chromosome, you will
    want to make this a small value.
}
  \item{pch}{
    The type of plotting symbol. See \code{\link{par}}.
}
  \item{colors}{
    A five-element character vector with the colors for: probes without
    change, probes that have a "gained" status, probes that have a
    "lost" status, the line that connects (smoothed values of) probes,
    the horizontal line at the 0 level.
  }

  \item{...}{Additional arguments; not used.}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{pSegment}}
}
\examples{

snowfallInit(universeSize = 4, typecluster = "SOCK")
fname <- list.files(path = system.file("data", package = "ADaCGH2"),
                     full.names = TRUE, pattern = "inputEx1")
tableChromArray <- inputDataToADaCGHData(filename = fname)

hs_mad.out <- pSegmentHaarSeg("cghData.RData",
                           "chromData.RData", merging = "MAD")

save(hs_mad.out, file = "hs_mad.out.RData", compress = FALSE)


pChromPlot(outRDataName = "hs_mad.out.RData",
                            cghRDataName = "cghData.RData",
                            chromRDataName = "chromData.RData",
                            posRDataName = "posData.RData",
                            probenamesRDataName = "probeNames.RData",
                            imgheight = 350)

}

\keyword{ hplot }
\keyword{ IO }
