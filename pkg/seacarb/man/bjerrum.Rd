\encoding{latin1}
\name{bjerrum}
\alias{bjerrum}
\title{Bjerrum plot}
\description{Plot the concentration of the various ionic forms of a molecule as a function of pH}
\usage{
bjerrum(K1=K1(), K2=NULL, K3=NULL, phmin=2, phmax=12, by=0.1, conc=1, 
     type="l", col="black",  ylab="Concentration (mol/kg)", add=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{K1}{First dissociation constant}
  \item{K2}{Second dissociation constant, default is NULL}
  \item{K3}{Third dissociation constant, default is NULL}
  \item{phmin}{Minimum pH value, default is 2}
  \item{phmax}{Maximum pH value, default is 12}
  \item{by}{Increment on the pH axis, default is 0.1}
  \item{conc}{concentration of molecule, default is 1}
  \item{type}{Type of plot, default is line}
  \item{col}{Color of plot, default is black}
  \item{ylab}{Label of Y axis, default is (mol/kg)}
  \item{add}{false:start new, true: add to current, default is false}
  \item{\dots}{Graphical parameters (see \code{\link{par}}) and any further arguments of plot, typically \code{\link{plot.default}}, may also be supplied as arguments to this function. Hence, the high-level graphics control arguments described under \code{\link{par}} and the arguments to \code{\link{title}} may be supplied to this function.}
}

%\value{XXX}

\references{
Zeebe, R. E. and Wolf-Gladrow D. A., 2001 \emph{CO2 in seawater: equilibrium, kinetics, isotopes}. Amsterdam: Elsevier, 346 pp.
}
             
\author{
Karline Soetaert \email{K.Soetaert@nioo.knaw.nl}
}

\seealso{
	\code{\link{matplot}}, \code{\link{par}}, \code{\link[seacarb]{speciation}}.
}

\details{
Note that the concentration is plotted in mol/kg only if conc is given is mol/kg
}

\examples{
## Plot the bjerrum plot for the carbonate system using the default values
bjerrum(K1(),K2(),main="DIC speciation",lwd=2) 
abline(v=-log10(K1()),col="grey")
mtext(side=3,at=-log10(K1()),"pK1")
abline(v=-log10(K2()),col="grey")
mtext(side=3,at=-log10(K2()),"pK2")
legend("left",lty=1:3,lwd=2,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for phosphate using the default values
bjerrum(K1p(),K2p(),K3p(),main="phosphate speciation",lwd=2)
legend("left",lty=1:4,lwd=2,legend=c(expression(H[3]~PO[4]),expression(H[2]~PO[4]^"-"),
expression(HPO[4]^"2-"),expression(PO[4]^"3-")))

## Plot the bjerrum plot for the carbonate system using the values other than the default ones,
## showing the effect of temperature
bjerrum(K1(T=25,S=35),K2(T=25,S=35),conc=1.3,main="effect of temperature" )
bjerrum(K1(T=0,S=35),K2(T=0,S=35),conc=1.3,add=TRUE,col="red")
legend("left",lty=1,col=c("black","red"),legend=c("T=25 oC","T=0 oC"))
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for the carbonate system using the values other than the default ones,
## showing the effect of salinity
bjerrum(K1(T=25,S=35),K2(T=25,S=35),conc=1.3,main="effect of salinity" )
bjerrum(K1(T=25,S=5),K2(T=25,S=5),conc=1.3,add=TRUE,col="blue")
legend("left",lty=1,col=c("black","blue"),legend=c("S=35","S=5"))
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))

## Plot the bjerrum plot for the carbonate system using the values other than the default ones,
## showing the effect of pressure
bjerrum(K1(P=0),K2(P=0),conc=1.3,main="effect of pressure" )
bjerrum(K1(P=300),K2(P=300),conc=1.3,add=TRUE,col="green")
legend("left",lty=1,col=c("black","green"),legend=c("P=0","P=300"),title="atm")
legend("right",lty=1:3,legend=c(expression(CO[2]),expression(HCO[3]^"-"),
	expression(CO[3]^"2-")))
}

\keyword{utilities}
