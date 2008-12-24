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

"K2p" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Phosphoric acid ---------------------
	#
	#   Guide to Best Practices in Ocean CO2 Measurements 2007 Chap 5 p 15
	#  (Dickson and Goyet): pH_T, mol/(kg-soln)
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnK2p <- -8814.715 / TK + 172.0883 - 27.927 * log(TK) +
		(-160.34 / TK + 1.3566) * sqrt(stp$S) + (0.37335 / TK - 0.05778) * stp$S
	res <- exp(lnK2p) * .KPcorr(stp, "K2p")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "K2p"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [HPO4--] / [H2PO4-]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(K2p()) # = -13.727