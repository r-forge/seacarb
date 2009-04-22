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

"K3p" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Phosphoric acid ---------------------
	#
	#   Guide to Best Practices in Ocean CO2 Measurements 2007 Chap 5 p 15
	#
	#   Ch.5 p. 15
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnK3p <- -3070.75 / TK - 18.141 + (17.27039 / TK + 2.81197) * sqrt(stp$S) +
		(-44.99486 / TK - 0.09984) * stp$S
	res <- exp(lnK3p) * .KPcorr(stp, "K3p")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "K3p"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [PO4---] / [HPO4--]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(K3p()) # = -20.24