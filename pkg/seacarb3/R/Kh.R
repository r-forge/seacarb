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

"Kh" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Kh (K Henry) ---------------------------------
	# Kh = [CO2*] / f(CO2) = [CO2] / p CO2   for CO2(g) <-> CO2(aq.)
	#
	# Weiss (1974)   [mol/kg/atm]
	# (Dickson, 1994 DOE chap 5 p 13)
	# f(CO2) is fugacity in atm. [ ] in mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)
	if (any(stp$P > 0)) warning("'Kh' calculated only for P = 0")

	TK <- .TK(stp$T)

	lnK0 <- 9345.17 / TK - 60.2409 + 23.3585 * log(TK / 100) +
		(0.023517 - 0.00023656 * TK + 0.0047036e-4 * TK^2) * stp$S
	res <- exp(lnK0)	# Correction for P > 0 ???

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kh"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[CO2*] / f(CO2)"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kh()) # = -3.5617