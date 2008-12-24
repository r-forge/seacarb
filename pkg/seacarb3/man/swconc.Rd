\encoding{latin1}
\name{swconc}
\alias{swconc}

\title{A seawater concentration object}

\description{An S3 object to represent the concentration of a chemical in seawater}

\usage{
swconc(x, substance = NULL, unit = getOption("seacarb.unit"),
	stp = NULL)
}
%- usage for print.swconc
\arguments{
  \item{x}{a vector of numerical data}
  \item{substance}{the name of the corresponding chemical substance}
  \item{unit}{unit of the concentration measurement}
  \item{stp}{an stp object representing salinity, temperature and pressure conditions}
}


\value{
  A 'swconc' object

}

\author{
Philippe Grosjean \email{phgrosjean@sciviews.org}
}


\examples{
# All conversion factors for standard STP conditions, without prefix
STP <- stp()
fromU <- rep(c("M", "m", "i"), each = 3)
toU <- rep(c("M", "m", "i"), times = 3)
convTable <- data.frame(to = rep(NA, 9), back = rep(NA, 9), check = rep(NA, 9))
for (i in 1:9) {
	convTable[i, 1] <- convert(swconc(1, unit = fromU[i], stp = STP), toU[i])
	convTable[i, 2] <- convert(swconc(1, unit = toU[i], stp = STP), fromU[i])
}
convTable$check <- convTable$to * convTable$back
row.names(convTable) <- paste(fromU, "->", toU)
convTable

# Same unit, but different prefixes
fromU <- rep(c("M", "mM", "uM", "nM"), each = 4)
toU <- rep(c("M", "mM", "uM", "nM"), times = 4)
convTable <- data.frame(to = rep(NA, 16), back = rep(NA, 16), check = rep(NA, 16))
for (i in 1:16) {
	convTable[i, 1] <- convert(swconc(1, unit = fromU[i], stp = STP), toU[i])
	convTable[i, 2] <- convert(swconc(1, unit = toU[i], stp = STP), fromU[i])
}
convTable$check <- convTable$to * convTable$back
row.names(convTable) <- paste(fromU, "->", toU)
convTable

# Conversions of concentrations according to changing STP conditions
STP1 <- stp(S = rep(c(35, 25), each = 6))
STP1
STP2 <- stp(S = rep(c(35, 25), each = 6),
	T = rep(c(5, 15, 25), times = 4),
	P = rep(c(0, 10, 0, 10), each = 3))
STP2
C1 <- swconc(rep(1, 12), unit = "molar", stp = STP1) 	# Initial conditions
C1

# Keep same unit (mol/L), but change T and/or P conditions
convert(C1, stp = STP2)

# Convert from molarity to molality AND for different T and/or P conditions
convert(C1, unit = "m", stp = STP2)

# It is an error to try changing salinity (not same seawater masses!)
convert(C1, unit = "m", stp = stp(S = 40))

# Example showing wrong and right way to change conc according to equilibrium!
# TODO
}

\keyword{utilities}
