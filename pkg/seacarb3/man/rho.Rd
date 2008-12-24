\encoding{latin1}
\name{rho}
\alias{rho}

\title{Density of seawater  (kg/m^3)}

\description{Calculates the density of seawater (\eqn{kg\ m^{-3}})}

\usage{
rho(x = NULL, S = 35, T = 25, P = 0)
}

\arguments{
  \item{x}{an stp object or an object containing an stp attribute, or \code{NULL}}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
}


\value{
  \item{rho}{Density of seawater (kg/m^3)}

}
\references{Millero, F. J. and Poisson, A. 1981 International one-atmosphere equation of state of seawater. \emph{Deep-Sea Research} \bold{28A}, 625-629.}

\author{
Aur\'elien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}


\examples{
rho(S = 35, T = 25, P = 0)
}

\keyword{utilities}
