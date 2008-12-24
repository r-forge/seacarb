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

"Kn" <- function (x = NULL, S = 35, T = 25, P = 0, ...)
{
	# Dissociation constant of ammonium on seawater scale - Millero 1995
	#
	# Kn = [H+] . [NH3] / [NH4+]
	# pH scale: originaly, seawater scale, but converted to total scale
	#

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnKn <- -6285.33 / TK + 0.0001635 * TK -0.25444 +
		(0.46532 - 123.7184 / TK) * sqrt(stp$S) +
		(-0.01992 + 3.17556 / TK) * stp$S
	res <- exp(lnKn) * .KPcorr(stp, "Kn") / Kconv(stp, ...)$Ktotal2sws

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kn"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [NH3] / [NH4+]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kn()) # = -21.3361???