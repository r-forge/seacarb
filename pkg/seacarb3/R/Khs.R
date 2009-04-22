# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
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

"Khs" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# First dissociation constant of H2S on total scale - Millero 1995
	#
	# Khs = [H+] . [HS-] / [H2S]
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnKhs <- 225.838 + 0.3449 * sqrt(stp$S) - 0.0274 * stp$S -
		13275.3 / TK - 34.6435 * log(TK)
	res <- exp(lnKhs) * .KPcorr(stp, "Khs")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Khs"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [HS-] / [H2S]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Khs()) # = -14.9908