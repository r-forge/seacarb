\encoding{latin1}
\name{Kspc}
\alias{Kspc}

\title{Solubility product of calcite (mol/kg)}

\description{Solubility product of calcite (mol/kg)}

\usage{
Kspc(x = NULL, S = 35, T = 25, P = 0)
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\value{
  \item{Kspc}{Solubility product of calcite (mol/kg)}

}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\references{Mucci, A. 1983 The solubility of calcite and aragonite in seawater at various salinities, temperature, and one atmosphere total pressure. \emph{American Journal of Science} \bold{283}: 780-799.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{Kspa}}.
}

\examples{
Kspc(S = 35, T = 25, P = 0)
log(Kspc()) # = -14.6659???
}
\keyword{utilities}
