\encoding{latin1}
\name{pHinsitu}
\alias{pHinsitu}

\title{pH at in situ temperature}

\description{pH in the sea with a different temperature that the one used to measure pH in the lab}

\usage{
pHinsitu(pH = 8.2, alk = 2.4e-3, Tinsitu = 20, Tlab = 25, S = 35, Pt = 0, Sit = 0)
}

\arguments{
  \item{pH}{pH measured in the laboratory}
  \item{alk}{total alkalinity (mol/kg-soln)}
  \item{Tinsitu}{In situ temperature in degrees Celsius}
  \item{Tlab}{Measurement temperature in degrees Celsius}
  \item{S}{Salinity}
  \item{Pt}{value of the concentration of total phosphate in mol/kg}
  \item{Sit}{the value of the total silicate in mol/kg}

}

\value{
  \item{pH}{pH at in situ temperature}
}

\references{
Hunter K. A., 1998. The temperature dependence of pH in surface seawater. \emph{Deep-Sea Research (Part I, Oceanographic Research Papers)} \bold{45}(11):1919-1930.
}

\author{
  Jean-Pierre Gattuso (\email{gattuso@obs-vlfr.fr})
}

%\seealso{	\code{\link{carb}}}

\examples{
# Trivial
pHinsitu(pH = 8.2, alk = 2.4e-3, Tinsitu = 25, Tlab = 25, S = 35)
# More interesting
pHinsitu(pH = 8.2, alk = 2.4e-3, Tinsitu = 10, Tlab = 25, S = 35)
}

\keyword{utilities}
