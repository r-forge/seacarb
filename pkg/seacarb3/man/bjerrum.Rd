\encoding{latin1}
\name{bjerrum}
\alias{bjerrum}

\title{Bjerrum plot}

\description{Plot the concentration of the various ionic forms of a molecule as a function of pH}

\usage{
bjerrum(conc = swconc(1, "ST", stp = stp()), pH = swpH(seq(0, 14, 0.1)),
	K1 = Ks(stp()), K2 = NULL, K3 = NULL, names = NULL, log = "", add = FALSE, ...)
}

\arguments{
  \item{conc}{A swconc object with concentrations to use for calculation}
  \item{pH}{A swPh object with pH to use for speciation}
  \item{K1}{First dissociation constant}
  \item{K2}{Second dissociation constant, default is NULL}
  \item{K3}{Third dissociation constant, default is NULL}
  \item{names}{names to use for the various ionic forms, e.g., \code{c("NH3", "NH4.")}.
    Please, follow naming convention and thus, use only letters, numbers, '.' or '_'.}
  \item{log}{plot concentration in log scale or not}
  \item{add}{do we add the plot to an existing graph?}
  \item{\dots}{further graphical arguments to pass to \code{plot}}
}

\value{
  Returns invisibly the concentrations of the various ionic forms for the
  different pHs. This function is mainly used for its side-effect of drawing a
  Bjerrum plot. Use \code{\link{speciation}} for just calculating the various
  ionic forms at given pHs.
}

\references{
  Zeebe, R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}

\author{
  Karline Soetaert (\email{K.Soetaert@nioo.knaw.nl}), slightly adapted by
  Philippe Grosjean (\email{Philippe.Grosjean@umh.ac.be})
}

\seealso{
  \code{\link{speciation}, \code{\link{matplot}}, \code{\link{par}}}
}

\examples{
# Bjerrum plot for the sulphate system (with default values)
bjerrum()
# see example(speciation) for further Bjerrum plots
}

\keyword{utilities}
