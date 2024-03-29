# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
"bor" <-
function(S = 35, T = 25, P = 0) {
	bor <- (416 * (S / 35)) * 1e-6   # (mol/kg), DOE94
	# PhG... and corrections for T and P???
	attr(bor, "unit") <- "mol/kg"
	return(bor)
}
