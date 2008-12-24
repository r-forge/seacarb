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
### TODO: pH conversion

"swpH" <- function (x, pHscale = getOption("seacarb.pHscale"), stp = NULL)
{
	if (inherits(x, "swpH")) 	# Convert pH scale, or calculate pHinsitu!
		return(convert(x, pHscale = pHscale, stp = stp))

	if (!is.null(stp))
		stp(x) <- stp 	# Set stp conditions for 'x'

	if (!is.numeric(x))
		stop("'x' must be numeric values")

	if (any(na.omit(x) < -1) || any(na.omit(x) > 15))
		stop("'x' must be -1 < pH < 15")

	attr(x, "pH scale") <- .checkpHscale(pHscale)

	attr(x, "class") <- c("swpH", "swobj")
	return(x)
}

"print.swpH" <- function (x, print.stp = TRUE, ...)
{
	if (!inherits(x, "swpH"))
		stop("'x' must be a 'swpH' object (seawater pH values)")

	X <- as.numeric(x)		# Strips out attributes and print values
	names(X) <- names(x)	# ... except names
	print(X, ...)

	cat("Seawater pH measured using ", attr(x, "pH scale"), "\n", sep = "")

	if (print.stp) {
		stp <- attr(x, "stp")
		if (!is.null(stp)) print(stp)
	}

	invisible(x)
}

"convert.swpH" <- function (x, pHscale = getOption("seacarb.pHscale"),
	stp = NULL, stp.recalc = TRUE, ...)
{
	if (!inherits(x, "swpH"))
		stop("'x' must be a 'swpH' object (seawater pH values)")

	pHscale <- .checkpHscale(pHscale)
	oldScale <- attr(x, "pH scale")

	if (!is.null(stp)) {
		oldStp <- attr(x, "stp")
		stp(x) <- stp

		# If stp was not set or stp.recalc == FALSE => set stp without recalc
		if (!is.null(oldStp) && stp.recalc) {
			# Possibly recalculate pH for new stp conditions
			### TODO...

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

			# Recalculate pH for new T and/or P conditions
			# (inspired by the original pHinsi() function)
			### TODO: requires also the alkalinity to allow the calc to be done!
		}
	}

	# It it a change in the scale?
	if (pHscale == oldScale) return(x)	# Nothing more to do!

	# Do the conversion TODO
	### TODO: get inspired by pHconv()

	# Update attributes
	attr(x, "pH scale") <- pHscale
	return(x)
}
