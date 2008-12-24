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
### TODO: assign to subset "[<-.swobj", especially another swobj
### Should be done object by object???

"[.swobj" <- function (x, i)
{
	# Allow subsetting 'swobj' objects without loosing attributes
	Names <- names(x)
	names(x) <- NULL
	Attr <- attributes(x)
	x <- unclass(x)[i]
	attributes(x) <- Attr
	names(x) <- Names[i]

	# stp must also be subsetted if it has more than 1 row
	stp <- attr(x, "stp")
	if (NROW(stp) > 1) stp(x) <- stp[i, ] else stp(x) <- stp

	return(x)
}

"as.data.frame.swobj" <- function (x, row.names = NULL, optional = FALSE, ...)
{
	Attr <- attributes(x)
	x <- unclass(x)
	res <- as.data.frame(x, row.names = row.names, optional = optional, ...)

	# Reinject attributes for each column in the resulting data.frame
	if (ncol(res) > 0)
		for (i in 1:ncol(res)) attributes(res[[i]]) <- Attr

	return(res)
}
