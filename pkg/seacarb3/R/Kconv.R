# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
# Slightly reworked by Ph. Grosjean, 2008, and renamed 'Kconv' for coherence
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

#--------------------------------------------------------------
# Conversion factors for converting dissociation constants
# from total pH scale to free pH scale (Ktotal2free)
# from total pH scale to seawater pH scale (Ktotal2sws)
# and from free pH scale to seawater scale (Kfree2sws)
# Kfree = Ktotal * Ktotal2free
# Ksws  = Ktotal * Ktotal2sws
# Ksws  = Kfree  * Kfree2sws
#--------------------------------------------------------------

### TODO: not same result as kconv() initial, but don't see the problem
# for Kconv(Kf.version = "dg")
"Kconv" <- function (x = NULL, S = 35, T = 25, P = 0, kf = "pf")
{
	stp <- stp(x, S, T, P)

	ST <- 0.14 / 96.062 / 1.80655 * stp$S    	# total sulfate
	Ks <- Ks(stp)
	FT <- 7e-5 * (stp$S / 35)                  	# total fluoride
	Kf <- Kf(stp, kf = kf)

	total2free <- 1 / (1 + ST / Ks)      		# Kfree = Ktotal * total2free

	Kf  = Kf*total2free       # convert Kf from total to free pH scale

	free2sws <- 1 + ST / Ks + FT / Kf      		# Ksws = Kfree * free2sws
	total2sws <- total2free * free2sws 			# Ksws = Ktotal * total2sws

	res <- list(
		Ktotal2sws = as.numeric(total2sws),
		Ktotal2free = as.numeric(total2free),
		Kfree2sws = as.numeric(free2sws))
	attr(res, "stp") <- stp
	attr(res, "Kf.version") <- attr(Kf, "version")

	return(res)
}
