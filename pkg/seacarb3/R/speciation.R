# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
# Refactored by Philippe Grosjean (2008) tu use new swobj objects
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA
#

"speciation" <-
	function (conc = swconc(1, "ST", stp = stp()), # total concentration
	pH = swpH(8),	   	            		       # pH used for speciation
	K1 = Ks(stp()), K2 = NULL, K3 = NULL, 	       # dissociation constants
	names = NULL, all.vars = FALSE, add.unit = FALSE)
{

	"univalent" <- function (conc, H, K1, Names)
	{
		# concentration for univalent species (e.g., NH3, NH4+)
		res <- list()
		C1 <- H / (K1 + H) # modified by HŽlo•se Lavigne (was C1=K1/(K1+H)
		L <- length(C1)
		if (L > length(conc)) conc[1:L] <- rep(conc, length.out = L)
		attr(conc, "substance") <- Names[[1]]
		res[[Names[1]]] <- conc * C1
		C2 <- 1 - C1
		attr(conc, "substance") <- Names[[2]]
		res[[Names[2]]] <- conc * C2
		return(res)
	}

	"bivalent" <- function (conc, H, K1, K2, Names)
	{
		# concentration for bivalent species (e.g., CO2, HCO3-, CO3--)
		res <- list()
		den <- H^2 + H * K1 + K1 * K2
		C1 <- H^2 / den
		L <- length(C1)
		if (L > length(conc)) conc[1:L] <- rep(conc, length.out = L)
		attr(conc, "substance") <- Names[[1]]
		res[[Names[1]]] <- conc * C1
		C2 <- H * K1 / den
		attr(conc, "substance") <- Names[[2]]
		res[[Names[2]]] <- conc * C2
		C3 <- 1 - C1 - C2
		attr(conc, "substance") <- Names[[3]]
		res[[Names[3]]] <- conc * C3
		return(res)
	}

	"trivalent" <- function (conc, H, K1, K2, K3, Names)
	{
		# concentration for trivalent species (H3PO4, H2PO4-, HPO4--, PO4---)
		res <- list()
		den <- H^3 + H^2 * K1 + H * K1 * K2 + K1 * K2 * K3
		C1 <- H^3 / den
		L <- length(C1)
		if (L > length(conc)) conc[1:L] <- rep(conc, length.out = L)
		attr(conc, "substance") <- Names[[1]]
		res[[Names[1]]] <- conc * C1
		C2 <- H^2 * K1 / den
		attr(conc, "substance") <- Names[[2]]
		res[[Names[2]]] <- conc * C2
		C3 <- H * K1 * K2 / den
		attr(conc, "substance") <- Names[[3]]
		res[[Names[3]]] <- conc * C3
		C4 <- 1 - C1 - C2 - C3
		attr(conc, "substance") <- Names[[4]]
		res[[Names[4]]] <- conc * C4
		return(res)
	}

	# Estimates the speciation of the various ionic forms of a molecule
	# as a function of pH

	# Check arguments (note: conc must provide STP conditions, and other STP
	# values must match this one, or be NULL)
	# Unit is given by conc too, and conversion possibly occurs for the rest
	if (!inherits(conc, "swconc")) conc <- swconc(conc, stp = stp())
	unit <- attr(conc, "unit")
	stp <- stp(conc)
	if (is.null(stp))
		stop("'stp' conditions must be defined for 'conc'")

	if (!inherits(pH, "swpH")) pH <- swpH(pH)
	if (!.compareStp(stp, stp(pH)))
		stop("'stp' conditions for 'conc' and 'pH' do not match")
	## TODO: check pH scale!! (=> otherwise, convert to total scale?)
	H <- as.numeric(10^-pH)   # proton concentration

	# Make sure names are character strings and we have (at least) four of them
	# also guess species names and add short unit
	if (is.null(names)) {
		if (add.unit) {
			Names <- .specNames(attr(conc, "substance"), unit)
		} else {
			Names <- .specNames(attr(conc, "substance"))
		}
	} else {
		l <- length(Names)
		if (l < 4) Names[(l+1):4] <- c("C1", "C2", "C3", "C4")[(l+1):4]
	}

	# Calculate speciation
	if (!inherits(K1, "swK"))
		stop("'K1' must be provided and must be a 'swK' object")
	if (!.compareStp(stp, stp(K1)))
		stop("'stp' conditions for 'conc' and 'K1' do not match")

	if (is.null(K2)) {
		res <- univalent(conc, H, as.numeric(K1), Names)
	} else if (is.null(K3)) {
		if (!inherits(K2, "swK"))
			stop("'K2' must be a 'swK' object")
		if (!.compareStp(stp, stp(K2)))
			stop("'stp' conditions for 'conc' and 'K2' do not match")
		res <- bivalent(conc, H, as.numeric(K1), as.numeric(K2), Names)
	} else {
		if (!inherits(K3, "swK"))
			stop("'K3' must be a 'swK' object")
		if (!.compareStp(stp, stp(K3)))
			stop("'stp' conditions for 'conc' and 'K3' do not match")
		res <- trivalent(conc, H, as.numeric(K1), as.numeric(K2),
			as.numeric(K3), Names)
	}

	# Create a 'speciation' object, inheriting from data.frame
	res <- as.data.frame(res)
	# Finalize the 'speciation' object
	attr(res, "unit") <- unit
	attr(res, "conc") <- conc
	attr(res, "pH") <- pH
	attr(res, "species") <- names(res)

	# Do we add also STP, pH and conc tot in the table?
	if (all.vars) {
		n <- nrow(res)
		Ctot <- attr(conc, "substance")
		if (is.null(Ctot)) Ctot <- "Total"
		if (add.unit)
			Ctot <- paste(Ctot, .shortUnit(attr(conc, "unit")), sep = "")
		if (length(conc) != n) conc[1:n] <- rep(conc, length.out = n)
		res[[Ctot]] <- conc
		if (length(pH) != n) pH[1:n] <- rep(pH, length.out = n)
		if (add.unit) {
			pHscale <- attr(pH, "pH scale")
			if (pHscale == "free scale") {
				res$pHF <- pH
			} else if (pHscale == "total scale") {
				res$pHT <- pH
			} else { # Seawater scale
				res$pHSWS <- pH
			}
			res$S <- stp$S
			res$T_C <- stp$T
			res$P_bar <- stp$P
		} else {
			res$pH <- pH
			res$S <- stp$S
			res$T <- stp$T
			res$P <- stp$P
		}
	}
	attr(res, "class") <- c("speciation", "data.frame")
	return(res)
}

"plot.speciation" <- function (x, y = NULL, type = "bjerrum", col = "black",
	lty = 1:4, lwd = 1, log = "", xlab = NULL, ylab = NULL, main = NULL,
	legend = NULL, pos = "right", add = FALSE, ...)
{
	if (!inherits(x, "speciation"))
		stop("'x' must be a 'speciation' object")

	Names <- names(x)

	# Depending on type, make a different graph, currently only "bjerrum"
	if (type == "bjerrum") {
		pH <- attr(x, "pH")
		species <- attr(x, "species")
		dat <- x[, species]

		if (is.null(xlab)) {
			pHscale <- attr(pH, "pH scale")
			if (pHscale == "free scale") {
				xlab <- expression(pH[F])
			} else if (pHscale == "total scale") {
				xlab <- expression(pH[T])
			} else { # Seawater scale
				xlab <- expression(pH[SWS])
			}
		}
		if (is.null(ylab)) {
			pUnit <- .prettyUnit(attr(x, "unit"))
			ylab <- bquote(paste("Concentration in ", .(pUnit)))
		}
		if (is.null(main)) {
			conc <- attr(x, "conc")

			# Check STP conditions
			stp <- attr(conc, "stp")
			if (NROW(stp) == 1) {
				cond <- bquote(paste(" (S = ", .(stp$S), ", T = ",
					.(stp$T) * degree, "C, P = ", .(stp$P), " bar)"))
			} else cond <- ""

			# Check substance
			subst <- attr(conc, "substance")
			if (is.null(subst)) {
				main <- bquote(paste("Bjerrum plot", .(cond)))
			} else {
				sbst <- .prettyLabel(subst)[[1]]
				main <- bquote(paste("Bjerrum plot for ", .(sbst), .(cond)))
			}
		}

		matplot(pH, dat, type = "l", col = col, lty = lty, lwd = lwd, log = log,
			xlab = xlab, ylab = ylab, main = main, add = add, ...)

		# Plot the legend
		if (!is.null(pos) && !add) {
			lgd <- .prettyLabel(species)
			legend(x = pos, inset = 0.05, legend = lgd, col = col, lty = lty,
				lwd = lwd, bg = "white")
		}
	}
}
