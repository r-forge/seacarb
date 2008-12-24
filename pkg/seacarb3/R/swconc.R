# Copyright (C) 2008, Philippe Grosjean
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

"swconc" <- function (x, substance = NULL, unit = getOption("seacarb.unit"),
	stp = NULL)
{
	if (!is.null(substance))	# Set/change substance for 'x'
		attr(x, "substance") <- as.character(substance[1])

	if (inherits(x, "swconc")) 	# Possibly convert unit!
		return(convert(x, unit = unit, stp = stp))

	if (!is.null(stp))
		stp(x) <- stp 	# Set stp conditions for 'x'

	if (!is.numeric(x))
		stop("'x' must be numeric values")

	attr(x, "unit") <- .checkUnit(unit)

	attr(x, "class") <- c("swconc", "swobj")
	return(x)
}

"print.swconc" <- function (x, print.stp = TRUE, ...)
{
	if (!inherits(x, "swconc"))
		stop("'x' must be a 'swconc' object (concentration)")

	X <- as.numeric(x)		# Strips out attributes and print values
	names(X) <- names(x)	# ... except names
	print(X, ...)

	substance <- attr(x, "substance")
	if (is.null(substance))
		txt <- "" else txt <- paste(substance, " in ", sep = "")

	unit <- attr(x, "unit")
	subunit <- sub("^.*mol", "mol", unit)
	if (subunit == "mol/kg-soln") {
		cat(txt, unit, " (molinity)\n", sep = "")
	} else if (subunit == "mol/kg-H2O") {
		cat(txt, unit, " (molality)\n", sep = "")
	} else if (subunit == "mol/L") {
		cat(txt, unit, " (molarity)\n", sep = "")
	} else stop ("Unknown unit: ", unit)

	if (print.stp) {
		stp <- attr(x, "stp")
		if (!is.null(stp)) print(stp)
	}

	invisible(x)
}

"convert.swconc" <- function (x, unit = getOption("seacarb.unit"),
	stp = NULL, stp.recalc = TRUE, ...)
{
	if (!inherits(x, "swconc"))
		stop("'x' must be a 'swconc' object (concentration)")

	unit <- .checkUnit(unit)
	oldUnit <- attr(x, "unit")

	if (!is.null(stp)) {
		oldStp <- attr(x, "stp")
		stp(x) <- stp

		# If stp was not set or stp.recalc == FALSE => set stp without recalc
		if (!is.null(oldStp) && stp.recalc) {
			# Possibly recalculate concentrations for new stp conditions
			# Warning! Use this feature only for conservative quantities like
			# total ion content (you MUST also recalculate constants Kxxx)!

			# Homogenize number of rows between the two objects
			n1 <- nrow(stp)
			n2 <- nrow(oldStp)
			n <- max(n1, n2)
			stp <- stp[rep(1:n1, length.out = n), ]
			oldStp <- oldStp[rep(1:n2, length.out = n), ]
			# Homogenize row.names
			row.names(stp) <- row.names(oldStp)

			# Salinity CANNOT be changed!
			if (!isTRUE(all.equal(stp$S, oldStp$S, scale = 1, tolerance = 0.1)))
				stop("Salinity cannot be changed")

			# Only molarities are affected by T and/or P changes
			if (sub("^.*/", "", oldUnit) == "L")
				x <- x * as.numeric(rho(stp) / rho(oldStp))
		}
	}

	# Is it a change in the unit?
	if (unit == oldUnit) return(x)	# Nothing more to do!

	# Do the conversion
	x <- x * .convUnit(oldUnit, unit, stp(x))

	# Update attributes
	attr(x, "unit") <- unit
	return(x)
}
