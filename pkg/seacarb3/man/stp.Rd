\encoding{latin1}
\name{stp}
\alias{stp}
\alias{stp<-}

\title{Salinity, temperature, pressure characteristics of seawater}

\description{The 'stp' object keeps seawater physical characteristics}

\usage{
stp(x = NULL, S = 35, T = 25, P = 0)
stp(x) <- value
}
%- usage for print.stp
\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, facultative}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{value}{new value to assign to x}
}


\value{
  \item{stp}{An object containing salinity, temperature, pressure
    characteristics of the medium (seawarter)}

}
\references{Millero, F. J. and Poisson, A. 1981 International one-atmosphere equation of state of seawater. \emph{Deep-Sea Research} \bold{28A}, 625-629.}

\author{
Philippe Grosjean \email{phgrosjean@sciviews.org}
}


\examples{
stp(S = 35, T = 25, P = 0)
}

\keyword{utilities}
