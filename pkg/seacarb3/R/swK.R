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
### TODO: convert.swK() using Kconv()

"print.swK" <- function (x, print.stp = TRUE, ...)
{
	if (!inherits(x, "swK"))
		stop("'x' must be a 'swK' object (seawater chemical equilibrium constant)")

	X <- as.numeric(x)		# Strips out attributes and print values
	names(X) <- names(x)	# ... except names
	print(X, ...)

	version <- attr(x, "version")
	if (is.null(version))
		version <- "" else version <- paste(", after ", version, sep = "")
	cat(attr(x, "K"), " (", attr(x, "unit"), ") = ", attr(x, "definition"),
		" (using pH ", attr(x, "pH scale"), version, ")\n", sep = "")

	if (print.stp) {
		stp <- attr(x, "stp")
		if (!is.null(stp)) print(stp)
	}

	invisible(x)
}
