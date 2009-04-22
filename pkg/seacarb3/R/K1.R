# Copyright (C) 2008 Jean-Pierre Gattuso and HŽlo•se Lavigne and Aurelien Proye
# Refactored by Ph. Grosjean, 2008
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

"K1" <- function (x = NULL, S = 35, T = 25, P = 0, k1k2 = "l")
{
	# --------------------- K1 with k1k2 = "l" ------------------------------
	# first acidity constant: K1 = [H+] . [HCO3-] / [CO2*]
	#
	# Mehrbach et al (1973) refit by Lueker et al. (2000).
	#
	# (Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
	# Dickson, Sabin and Christian , 2007, Chapter 5, p. 13)
	#
	# pH-scale: 'total'. mol/kg-soln

	# --------------------- K1 with k1k2 = "r" ------------------------------
	#
	# (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 14)
	#   pH-scale: 'total'. mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)
	if (!identical(k1k2, "l") && !identical(k1k2, "r"))
		stop("'k1k2' must be 'l' for Lueker et al., 2000 or 'r' for Roy et al., 1993")

	TK <- .TK(stp$T)

	if (k1k2 == "l") {
		logK1lue <- -3633.86 / TK + 61.2172 - 9.67770 * log(TK) +
			0.011555 * stp$S - 0.0001152* stp$S^2
		res <- 10^logK1lue * .KPcorr(stp, "K1")
		version <- "Lueker et al., 2000"
	} else {
		lnK1roy <- 2.83655 - 2307.1266 / TK - 1.5529413 * log(TK) -
			(0.20760841 + 4.0484 / TK) * sqrt(stp$S) + 0.08468345 * stp$S -
			0.00654208 * stp$S^(3/2) + log(1 - 0.001005 * stp$S)
		res <- exp(lnK1roy) * .KPcorr(stp, "K1")
		version <- "Roy et al., 1993"
	}

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "K1"	# Name of the constant
	attr(res, "k1k2") <- k1k2
	attr(res, "version") <- version
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"

	attr(res, "definition") <- "[H+] . [HCO3-] / [CO2*]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(K1(k1k2 = "r")) # = -13.4847