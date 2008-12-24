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

"Ksi" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# Dissociation constant of Si(OH)4 on total scale
	# DOE 1994 - chap 5 p 17
	#
	# Ksi = [H+] . [SiO(OH)3-] / [Si(OH)4]

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)
	iom0 <- .iom0(stp$S)

	lnKsi <- -8904.2 / TK + 117.385 - 19.334 * log(TK) +
		(-458.79 / TK + 3.5913) * sqrt(iom0) + (188.74 / TK - 1.5998) * iom0 +
		(-12.1652 / TK + 0.07871) * iom0^2 + log(1 - 0.001005 * stp$S)
	res <- exp(lnKsi) * .KPcorr(stp, "Ksi")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Ksi"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [SiO(OH)3-] / [Si(OH)4]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Ksi()) # = -21.61