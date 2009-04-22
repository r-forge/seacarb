# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
# Refactoring by Ph. Grosjean, 2008
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

"bjerrum" <-
function (conc = swconc(1, "ST", stp = stp()), # total concentration (molality)
	pH = swpH(seq(0, 14, 0.1)),	   	           # pH used for speciation
	K1 = Ks(stp()), K2 = NULL, K3 = NULL, 	   # dissociation constants
	names = NULL, log = "", add = FALSE, ...)
{
	# Creates a bjerrum plot (create a speciation object and plot it)
	# Use a better default for pH than speciation() alone
	sp <- speciation(conc = conc, pH = pH, K1 = K1, K2 = K2, K3 = K3,
		names = names, all.vars = TRUE)
	plot(sp, log = log, add = add, ...)
	return(invisible(sp))
}

# bjerrum() simply for the default Bjerrum plot (sulfate)
