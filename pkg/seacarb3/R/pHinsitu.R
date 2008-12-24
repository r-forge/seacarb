# Copyright (C) 2003 Jean-Pierre Gattuso and Aurelien Proye
# Slightly reworked and renamed 'pHinsitu' by Ph. Grosjean, 2008
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
"pHinsitu" <-
function (pH = 8.2, alk = 2.4e-3, Tinsitu = 20, Tlab = 25, S = 35, Pt = 0, Sit = 0)
{
	# according to Hunter (1998)
	#si TA = 0 alors calculer TA = 660+47.6 * S
	dat1 <- carb(flag = 8, var1 = pH, var2 = alk, S = S, T = Tlab,
		P = 0, Pt = Pt, Sit = 0)
	dat2 <- carb(flag = 15, var1 = alk, var2 = dat1$DIC, S = S,
		T = Tinsitu, P = 0, Pt = Pt, Sit = Sit)
	### TODO PhG: add the various required attributes
	return(dat2$pH)
}
