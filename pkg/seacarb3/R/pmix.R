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
"pmix" <-
function(flag, var1, var2, pCO2s, wf, S=35, T=20, P=0, Pt=0, Sit=0){
	sys <- 0
	if (sys==0) {
		ws <- 1 * wf
		wi <- 1 - ws
		ci <- carb(flag=flag, var1=var1, var2=var2, S=S ,T=T, P=P, Pt=Pt, Sit=Sit)
		cs <- carb(24, var1=pCO2s, var2=ci$ALK, S=S ,T=T, P=P, Pt=Pt, Sit=Sit) # carbonate chemistry of high-CO2 component
		alkf <- ci$ALK  # final alkalinity
		dicf <- (ci$DIC*wi + cs$DIC*ws)/(wi+ws)	# final dic
		cf <- carb(flag=15, var1=alkf, var2=dicf, S=S ,T=T, P=P, Pt=Pt, Sit=Sit)	 # final carbonate chemistry
		co <- as.data.frame(c("pmix-closed-initial", rep("pmix-closed-final", nrow(cf)))) # comment column
		out <- rbind(ci, cf)
		out <- cbind(co, out)
		names(out)[1] <- "comment"
		return(out[,1:15])
	}
	else { # sys=1 is useless as the final (after equilibration) carbonate chemistry is identical to the initial one
		return("invalid parameter(s)")
	}
}
