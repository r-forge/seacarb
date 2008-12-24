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
# Conversion factors for converting pH
# from total pH scale to free pH scale (pHtotal2free)
# from total pH scale to seawater pH scale (pHtotal2sws)
# and from free pH scale to seawater scale (pHfree2Ssws)
# pHfree = pHtotal + pHtotal2free
# pHsws  = pHtotal + pHtotal2SWS
# pHsws  = pHfree  + pHfree2SWS
#--------------------------------------------------------------

"pHconv" <- function (x = NULL, S = 35, T = 25, P = 0, kf = "pf")
{
	cc <- Kconv(x = x, S = S, T = T, P = P, kf = kf)

	res <- list(
		pHtotal2sws =  log10(1 / cc$Ktotal2sws),
		pHtotal2free = log10(1 / cc$Ktotal2free),
        pHfree2sws =   log10(1 / cc$Kfree2sws))

	attr(res, "stp") <- attr(cc, "stp")
	attr(res, "kf") <- attr(cc, "kf")

	return(res)
}
