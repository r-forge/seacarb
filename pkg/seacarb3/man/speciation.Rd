\encoding{latin1}
\name{speciation}
\alias{speciation}

\title{ionic forms as a function of pH}

\description{Estimate the concentration of the various ionic forms of a molecule as a function of pH}

\usage{
% also usage for plot.speciation
speciation(conc = swconc(1, "ST", stp = stp()), pH = swpH(8),
	K1 = Ks(stp()), K2 = NULL, K3 = NULL,
	names = NULL, all.vars = FALSE, add.unit = FALSE)
}

\arguments{
  \item{conc}{A swconc object with concentrations to use for calculation}
  \item{pH}{A swPh object with pH to use for speciation}
  \item{K1}{First dissociation constant}
  \item{K2}{Second dissociation constant, default is NULL}
  \item{K3}{Third dissociation constant, default is NULL}
  \item{names}{names to use for the various ionic forms, e.g., \code{c("NH3", "NH4.")}.
    Please, follow naming convention and thus, use only letters, numbers, '.' or '_'.}
  \item{all.vars}{do we also add stp, pH and conc tot in the table?}
  \item{add.unit}{do we add units to the names of the columns in the table}
}

\value{The function returns a table containing the following items:
  % TODO: update this!
  \item{C1}{ionic form 1, univalent, bivalent and trivalent molecules}
  \item{C2}{ionic form 2, univalent, bivalent and trivalent molecules}
  \item{C3}{ionic form 3, bivalent and trivalent molecules}
  \item{C4}{ionic form 4, trivalent molecules}
}

\references{
  Zeebe, R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}

\author{
 Karline Soetaert (\email{K.Soetaert@nioo.knaw.nl}), slightly adapted by
 Philippe Grosjean (\email{Philippe.Grosjean@umh.ac.be})
}

\seealso{
  \code{\link{bjerrum}}
}

\examples{
# Sulfate (univalent, using free pH scale)
SpS <- speciation(swconc(0.14/96.062/1.80655*35*1000, "ST", unit = "mmol/kg-soln",
	stp = stp()), pH = swpH(seq(0, 14, 0.1), "F"), K1 = Ks(stp()))
head(SpS)
plot(SpS, log = "y")	# Log axis for Y
plot(SpS)
# superpose curves for T = 10
SpS10 <- speciation(swconc(0.14/96.062/1.80655*35*1000, "ST", unit = "mmol/kg-soln",
	stp = stp(T = 10)), pH = swpH(seq(0, 14, 0.1), "F"), K1 = Ks(stp(T = 10)))
plot(SpS10, col = "red", add = TRUE)
legend("topright", inset = 0.05, col = c("red", "black"), lty = 1, bg = "white",
	legend = c(expression(10*degree*C), expression(25*degree*C)))

# Carbonate (bivalent ion)
SpC <- speciation(swconc(2.04, "CT", unit = "mmol/kg-soln", stp = stp()),
	pH = swpH(seq(0, 14, 0.1), "T"), K1 = K1(stp()), K2 = K2(stp()))
head(SpC)
plot(SpC, log = "y")
plot(SpC)
# superpose curves using Roy et al., 1993 algorithm
SpCR <- speciation(swconc(2.04, "CT", unit = "mmol/kg-soln", stp = stp()),
	pH = swpH(seq(0, 14, 0.1), "T"), K1 = K1(stp(), version = "r"),
	K2 = K2(stp(), version = "r"))
plot(SpCR, col = "red", add = TRUE)
legend("left", inset = 0.05, col = c("red", "black"), lty = 1, bg = "white",
	legend = c("Roy et al., 1993", "Lueker et al., 2000"))

# Phosphate (trivalent ion)
SpP <- speciation(swconc(1, "PT", unit = "umol/kg-soln", stp = stp()),
	pH = swpH(seq(0, 14, 0.1), "T"), K1 = K1p(stp()), K2 = K2p(stp()), K3 = K3p(stp()))
head(SpP)
plot(SpP)

# Fluorure
# TODO

# Borate
# TODO

# Ammonium
# TODO

# Silicate
# TODO
}

\keyword{utilities}
