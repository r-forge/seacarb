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

"Kb" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	#  Kb = [H+] . [B(OH)4-] / [B(OH)3]
	#
	#  (Dickson, 1990 in Guide to Best Practices in Ocean CO2 Measurements 2007)
	#  pH-scale: 'total'. mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnKb <- (-8966.90 - 2890.53 * sqrt(stp$S) - 77.942 * stp$S +
		1.728 * stp$S^(3/2) - 0.0996 * stp$S^2) / TK + 148.0248 +
		137.1942 * sqrt(stp$S) + 1.62142 * stp$S +
		(-24.4344 - 25.085 * sqrt(stp$S) - 0.2474 * stp$S) * log(TK) +
		0.053105 * sqrt(stp$S) * TK
	res <- exp(lnKb) * .KPcorr(stp, "Kb")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kb"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [B(OH)4-] / [B(OH)3]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kb()) # = -19.7964