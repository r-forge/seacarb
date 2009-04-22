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

"Ks" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# Dickson (1990) Equilibrium constant for HSO4- <=> H+ + SO4--
	#
	# Ks = [H+]free [SO4--] / [HSO4-]
	# pH scale: free scale
	#
	# the term log(1 - 0.001005 * S) converts from mol/kg-H2O to mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)
	iom0 <- .iom0(stp$S)

	lnKs <- -4276.1 / TK + 141.328 -23.093 * log(TK) +
		(-13856 / TK + 324.57 - 47.986 * log(TK)) * sqrt(iom0) +
		(35474 / TK - 771.54 + 114.723 * log(TK)) * iom0 +
		-2698 / TK * iom0^(3/2) + 1776 / TK * iom0^2 +
	    log(1 - 0.001005 * stp$S)
	res <- exp(lnKs) * .KPcorr(stp, "Ks")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Ks"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "free scale"
	attr(res, "definition") <- "[H+] . [SO4--] / [HSO4-]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Ks()) # = -2.30