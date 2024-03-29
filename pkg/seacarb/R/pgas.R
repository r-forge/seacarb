# Copyright (C) 2008 Jean-Pierre Gattuso
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
"pgas" <- function(flag, var1, var2, pCO2g, S=35, T=20, P=0, Pt=0, Sit=0, k1k2='l', kf='pf', pHscale="T"){
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, pHscale=pHscale)
		alkf <- ci$ALK
		cf <- carb(flag=24,var1=pCO2g, var2=alkf, S=S, T=T, P=P,  Pt=Pt, Sit=Sit, k1k2=k1k2, kf=kf, pHscale=pHscale)		
		co <- as.data.frame(c("pgas-initial", rep("pgas-final", nrow(cf))))
		out <- rbind(ci, cf)
		out <- cbind(co, out)
		names(out)[1] <- "comment"
		return(out)
}