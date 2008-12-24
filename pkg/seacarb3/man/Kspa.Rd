\encoding{latin1}
\name{Kspa}
\alias{Kspa}

\title{Solubility product of aragonite (mol/kg)}

\description{Solubility product of aragonite (mol/kg)}

\usage{
Kspa(x = NULL, S = 35, T = 25, P = 0)
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{Kspa}{Solubility product of aragonite (mol/kg)}

}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\references{Mucci, A. 1983 The solubility of calcite and aragonite in seawater at various salinities, temperature, and one atmosphere total pressure. \emph{American Journal of Science} \bold{283}: 780-799.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{Kspc}}.
}

\examples{
Kspa(S = 35, T = 25, P = 0)
log(Kspa()) # = -14.2491???
}
\keyword{utilities}
