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

"Kf" <- function (x = NULL, S = 35, T = 25, P = 0, kf = "pf")
{
	#---------------------- Kf Perez and Fraga ---------------------------
	#  Kf = [H+] . [F-] / [HF]
	#
	# Perez and Fraga, 1987 in Guide to the Best Practices for Ocean CO2 Measurements
	# Dickson, Sabin and Christian , 2007, Chapter 5, p. 14)
	#
	# --------------------- Kf Dickson and Goyet -------------------------
	# (Dickson and Riley, 1979 in Dickson and Goyet,
	# 1994, Chapter 5, p. 14)
	#
	# pH scale: 'total', mol/kg-soln
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)
		if (!identical(kf, "pf") && !identical(kf, "dg"))
		stop("'kf' must be 'pf' for Perez and Fraga, 1987 or 'dg' for Dickson and Goyet, 1979")

	TK <- .TK(stp$T)

	if (kf == "pf") {
		lnKf <- 874 / TK - 9.68 + 0.111 * sqrt(stp$S)
		version <- "Perez and Fraga, 1987"
	} else {
		iom0 <- .iom0(stp$S)
		ST <- 0.14 / 96.062 / 1.80655 * stp$S   		 # total sulfate
		lnKf <- 1590.2 / TK - 12.641 + 1.525 * sqrt(iom0) +
			log(1 - 0.001005 * stp$S) + log(1 + ST / Ks(stp))
		version <- "Dickson and Goyet, 1979"
	}
	res <- exp(lnKf) * .KPcorr(stp, "Kf")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kf"	# Name of the constant
	attr(res, "kf") <- kf
	attr(res, "version") <- version
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [F-] / [HF]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kf(kf = "dg")) # = -5.80