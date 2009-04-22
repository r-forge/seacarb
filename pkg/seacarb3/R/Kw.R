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

"Kw" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Kwater -------------------------------------
	#
	# Millero (1995)(in Guide to Best Practices in Ocean CO2 Measurements (2007, Chapter 5, p.16))
	# Kw in mol/kg-soln. pH scale: total scale.
	# Kw = [H+] . [OH-]


	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	lnKw <- -13847.26 / TK + 148.9652 - 23.6521 * log(TK) +
		(118.67 / TK - 5.977 + 1.0495 * log(TK)) * sqrt(stp$S) - 0.01615 * stp$S
	res <- exp(lnKw) * .KPcorr(stp, "Kw")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kw"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[H+] . [OH-]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kw()) # = -30.434