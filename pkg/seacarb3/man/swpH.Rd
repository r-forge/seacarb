\encoding{latin1}
\name{swpH}
\alias{swpH}

\title{A seawater pH object}

\description{An S3 object to represent the pH in seawater}

\usage{
swpH(x, pHscale = getOption("seacarb.pHscale"), stp = NULL)
}
%- usage for print.swpH and convert.swpH
\arguments{
  \item{x}{a vector of numerical data}
  \item{pHscale}{the scale to use for pH: free scale, total scale or sewater scale}
  \item{stp}{an stp object representing salinity, temperature and pressure conditions}
}


\value{
  A 'swpH' object

}

\author{
Philippe Grosjean \email{phgrosjean@sciviews.org}
}


\examples{
# TODO
}

\keyword{utilities}
