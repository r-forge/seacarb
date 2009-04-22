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

"K2" <- function (x = NULL, S = 35, T = 25, P = 0, k1k2 = "l")
{
	# --------------------- K2 with k1k2 = "l" ------------------------------
	# second acidity constant: K2 = [H+] . [CO3--] / [HCO3-]
	#
	# Mehrbach et al (1973) refit by Lueker et al. (2000).
	#
	# (Lueker  et al., 2000 in Guide to the Best Practices for Ocean CO2 Measurements
	# Dickson, Sabin and Christian , 2007, Chapter 5, p. 14)
	#
	# pH-scale: 'total'. mol/kg-soln

	# --------------------- K2 with k1k2 = "r" ------------------------------
	#
	# (Roy et al., 1993 in Dickson and Goyet, 1994, Chapter 5, p. 15)
	#   pH-scale: 'total'. mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)
	if (!identical(k1k2, "l") && !identical(k1k2, "r"))
		stop("'k1k2' must be 'l' for Lueker et al., 2000 or 'r' for Roy et al., 1993")

	TK <- .TK(stp$T)

	if (k1k2 == "l") {
		logK2lue <- -471.78 / TK - 25.9290 + 3.16967 * log(TK) +
			0.01781 * stp$S - 0.0001122 * stp$S^2
		res <- 10^logK2lue * .KPcorr(stp, "K2")
		version <- "Lueker et al., 2000"
	} else {
		lnK2roy <- -9.226508 - 3351.6106 / TK - 0.2005743 * log(TK) +
			(-0.106901773 - 23.9722 / TK) * sqrt(stp$S) +
			0.1130822 * stp$S - 0.00846934 * stp$S^(3/2) +
			log(1 - 0.001005 * stp$S)
		res <- exp(lnK2roy) * .KPcorr(stp, "K2")
		version <- "Roy et al., 1993"
	}

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "K2"	# Name of the constant
	attr(res, "k1k2") <- k1k2
	attr(res, "version") <- version
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"

	attr(res, "definition") <- "[H+] . [CO3--] / [HCO3-]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(K1(k1k2 = "r")) # = -20.5504