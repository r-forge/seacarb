\encoding{latin1}
\name{Kw}
\alias{Kw}

\title{Ion product of water (mol^2/kg^2)}

\description{Ion product of water (mol^2/kg^2)}

\usage{
Kw(x = NULL, S = 35, T = 25, P = 0)
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{Kw}{Ion product of water (mol2/kg2), pHscale = totale scale}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\references{Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica et Cosmochimica Acta} \bold{59} 661-677.}

\author{
Héloïse Lavigne, Aurélien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\examples{
Kw(S = 35, T = 25, P = 0)
log(Kw()) # = -30.434
}

\keyword{utilities}
