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

"Kspc" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# --------------------- Kspa (calcite) ----------------------------
	#
	# apparent solubility product of aragonite
	#
	#  Kspa = [Ca++]T [CO3--]T
	#
	#  where []T refers to the equilibrium total
	# (free + complexed) ion concentration.
	#
	#  Mucci 1983 mol/kg-soln

	# Check and rework arguments
	stp <- stp(x, S, T, P)

	TK <- .TK(stp$T)

	logKspc <- -171.9065 - 0.077993 * TK + 2839.319 / TK + 71.595 * log10(TK) +
		(-0.77712 + 0.0028426 * TK + 178.34 / TK) * sqrt(stp$S) -
		0.07711 * stp$S + 0.0041249 * stp$S^(3/2)
	res <- 10^(logKspc) * .KPcorr(stp, "Kspc")

	# Construct the 'swK' object together with its metadata
	attr(res, "K") <- "Kspc"	# Name of the constant
	attr(res, "unit") <- "mol/kg-soln"
	attr(res, "pH scale") <- "total scale"
	attr(res, "definition") <- "[Ca++] . [CO3--]"
	attr(res, "stp") <- stp
	attr(res, "class") <- c("swK", "swobj")
	return(res)
}

# log(Kspc()) # = -14.6659???