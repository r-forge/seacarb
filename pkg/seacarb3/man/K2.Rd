\encoding{latin1}
\name{K2}
\alias{K2}

\title{Second dissociation constant of carbonic acid (mol/kg)}

\description{Second dissociation constant of carbonic acid (mol/kg)}

\usage{
K2(x = NULL, S = 35, T = 25, P = 0, k1k2 = "l")
}

\arguments{
  \item{x}{an stp object, or an object containing an stp attribute, or NULL}
  \item{S}{Salinity, default is 35}
  \item{T}{Temperature in degrees Celsius, default is 25oC}
  \item{P}{Hydrostatic pressure in bar (surface = 0), default is 0}
  \item{k1k2}{"l" for using K1 and K2 from Lueker et al. and "r" for using K1 and K2 from Roy and al. , default is "l"}
 }

\value{
  \item{K2}{Second dissociation constant of carbonic acid (mol/kg),  pH scale = total scale}
}

\details{Note that the arguments can be given as a unique number or as vectors. If the lengths of the vectors are different, the longer vector is retained and only the first value of the other vectors is used. It can therefore be critical to use vectors of the same length.}

\references{The Lueker et al. (2000) constant is recommended by Guide to Best Practices for Ocean CO2 Measurements (2007). The Roy et al. (1993) constants is recommended by DOE (1994).

Dickson A. G., Sabine C. L. and Christian J. R., 2007 Guide to best practices for ocean CO2 measurements. \emph{PICES Special Publication} \bold{3}, 1-191.

DOE 1994 \emph{Handbook of methods for the analysis of the various parameters of the carbon dioxide system in sea water}. ORNL/CDIAC-74. Oak Ridge,Tenn.: Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory.

Lueker, T. J. Dickson, A. G. and Keeling, C. D. 2000 Ocean pCO2 calculated from dissolved inorganic carbon, alkalinity, and equations for K1 and K2: validation based on laboratory measurements of CO2 in gas and seawater at equilibrium. \emph{Marine Chemistry} \bold{70} 105-119.

Roy, R. N. Roy, L. N. Vogel, K. M. Porter-Moore, C. Pearson, T. Good C. E., Millero F. J. and Campbell D. M., 1993. The dissociation constants of carbonic acid in seawater at salinities 5 to 45 and temperatures 0 to 45�C. \emph{Marine Chemistry} \bold{44}, 249-267.
}

\author{
H�lo�se Lavigne, Aur�lien Proye and Jean-Pierre Gattuso \email{gattuso@obs-vlfr.fr}
}

\seealso{
	\code{\link{K1}}.
}

\examples{
K2(S = 35, T = 25, P = 0)
log(K1(k1k2 = "r")) # = -20.5504
}

\keyword{utilities}
