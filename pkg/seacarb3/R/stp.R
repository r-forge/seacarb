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

"stp" <- function (x = NULL, S = 35, T = 25, P = 0)
{
	# If 'x' is an stp object, return it immediately
	if (inherits(x, 'stp')) return(x)

	# Note: stp is also used to extract stp attribute from an object!
	stp <- attr(x, "stp")
	if (!is.null(stp)) return(stp)

	# If 'x' is a data.frame or a list, create the object from 'S', 'T' and 'P' components
	if (inherits(x, c("list", "data.frame"))) {
		if (!is.null(x$T)) T <- as.numeric(x$T)
		if (!is.null(x$P)) P <- as.numeric(x$P)
		if (!is.null(x$S)) S <- as.numeric(x$S)
	} else if (!is.null(x)) return(NULL)	# 'x' provided, but no stp data

	# Note: limits are arbitrarily set, still must be checked!
	if (!is.vector(S) || !is.numeric(S))
	    stop("'S' must be a vector of numeric values")
	if (any(is.na(S)) || any(S < 0) || any(S > 50))
	    stop("'S' cannot have missing values, or 0 < values < 50")
    if (!is.vector(T) || !is.numeric(T))
	    stop("'T' must be a vector of numeric values")
	if (any(is.na(T)) || any(T < 0) || any(T > 1000))
	    stop("'T' cannot have missing or negative values, or values higher than 1000 degrees C")
    if (!is.vector(P) || !is.numeric(P))
	    stop("'P' must be a vector of numeric values")
	if (any(is.na(P)) || any(P < 0) || any(P > 1000))
	    stop("'P' cannot have missing or negative values, or values higher than 1000 bar")

	res <- data.frame(S = S, T = T, P = P)
	attr(res,"class") <- c("stp", "data.frame")
	return(res)
}

"stp<-" <- function (x, value)
{
	# Special case: assigning NULL 'value' erases the 'stp' argument
	if (is.null(value)) {
		attr(x, 'stp') <- NULL
		return(x)
	}

	# Check that 'value' is an 'stp' object
	if (!inherits(value, "stp"))
		stop("You can assign a 'stp' object only")

	# Does value contain more than one row?
	if (NROW(value) > 1) {
		# Check if value is of same length as the object
		if (nrow(value) != NROW(x))
			stop("Can only assign a 'stp' object with 1, or the same number of rows as the target object")

		# If all S, T, and P values are identical, simplify to a one-row object
		"isConstant" <- function (x)
		{
			if (any(is.na(x))) return(FALSE)
			Rng <- range(x)
			return(Rng[1] == Rng[2])
		}
		if (isConstant(value$S) && isConstant(value$T) && isConstant(value$P)) {
			value <- value[1, ]
			row.names(value) <- NULL
		} else {
			# Make sure to adapt row.names to names of the object items in 'x'
			row.names(value) <- names(x)
		}
	} else row.names(value) <- NULL 	# Only one row, do not name it

	# Assign 'stp' as attribute to 'x'
	attr(x, "stp") <- value
	return(x)
}

"print.stp" <- function (x, ...)
{
	if (!inherits(x, "stp"))
		stop("'x' must be a 'stp' object (Salinity, Temperature, Pressure definition)")

	# If x contain unique values, use a condensed presentation
	if (length(x$S) == 1 && length(x$T) == 1 && length(x$P) == 1) {
		cat("Seawater S = ", x$S, ", T = ", x$T, ", P = ", x$P, "\n", sep = "")
	} else {	# Create a table of values
		cat("Seawater characteristics:\n")
		stp <- x
		class(stp) <- "data.frame"
		print(stp)
	}
	invisible(x)
}
