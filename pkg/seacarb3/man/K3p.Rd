\encoding{latin1}
\name{K3p}
\alias{K3p}

\title{Third dissociation constant of phosphoric acid (mol/kg)}

\description{Third dissociation constant of phosphoric acid (mol/kg)}

\usage{
K3p(x = NULL, S = 35, T = 25, P = 0)
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\value{
 \item{K3p}{Third dissociation constant of phosphoric acid (mol/kg),,pHscale = totale scale}
}

\references{
Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.}

\author{
H�lo�se Lavigne, Aur�lien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K1p}}, \code{\link{K2p}}.
}

\examples{
K3p(S = 35, T = 25, P = 0)
log(K3p()) # = -20.24
}

\keyword{utilities}
