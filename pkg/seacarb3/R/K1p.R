# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
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

"K1p" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Phosphoric acid ---------------------
	#
	#   (DOE, 1994)  (Dickson and Goyet): pH_T, mol/(kg-soln)
	#   Ch.5 p. 16
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnK1p <- -4576.752 / TK + 115.525 - 18.453 * log(TK) +
		(-106.736 / TK + 0.69171) * sqrt(stp$S) + (-0.65643 / TK - 0.01844) * stp$S
	res <- exp(lnK1p) * .KPcorr(stp, "K1p")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "K1p"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [H2PO4-] / [H3PO4]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

#log(K1p()) # = -3.71