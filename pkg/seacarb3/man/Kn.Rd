\encoding{latin1}
\name{Kn}
\alias{Kn}

\title{Dissociation constant of ammonium (mol/kg)}

\description{Dissociation constant of ammonium on the total scale (mol/kg)}

\usage{
Kn(x = NULL, S = 35, T = 25, P = 0, ...)
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{\dots}{Further arguments passed to the function}
}

\value{
	\item{Kn}{Dissociation constant of ammonium on the total scale (mol/kg)}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\references{
Millero F. J., 1995. Thermodynamics of the carbon dioxide system in the oceans. \emph{Geochimica Cosmochimica Acta} 59: 661-677.
}

\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

\examples{
Kn(S = 35, T = 25, P = 0)
log(Kn()) # = -21.3361???
}

\keyword{utilities}
