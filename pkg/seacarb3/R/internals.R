# Copyright (C) 2008 Philippe Grosjean
# Adapted from code by Jean-Pierre Gattuso and Aurelien Proye, 2003
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

# Internal functions used by seacarb
# Created by Ph. Grosjean, 2008, after refactoring original seacarb code

".TK" <- function (T)
{
	# Calculate TK in degrees Kelvin, from T in degrees Celsius
	return(T + 273.15)
}

".iom0" <- function (S)
{
	# Calculate ionic strength of seawater from salinity S
	# according to Dickson 1990 in DOE 1994 (see ?Ks for instance)
	return(19.924 * S / (1000 - 1.005 * S))
}

".KPcorr" <- function (stp, K = "K1")
{
	# Ksi, Kn and Khs added by Ph. Grosjean 2008 from code by Karline Soetaert
	# Pressure Correction from Millero 1995
	# Typos corrections from Lewis and Wallace (CO2SYS)
	# Calculation of pression correction for K with K being:
	# K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi, Kn or Khs
	# i is the corresponding vector index in the tables
	i <- switch(K,
		K1 = 1,
		K2 = 2,
		Kb = 3,
		Kw = 4,
		Ks = 5,
		Kf = 6,
		Kspc = 7,
		Kspa = 8,
		K1p = 9,
		K2p = 10,
		K3p = 11,
		Ksi = 12,
		Kn = 13,
		Khs = 14,
		stop("'K' must be K1, K2, Kb, Kw, Ks, Kf, Kspc, Kspa, K1p, K2p, K3p, Ksi, Kn or Khs"))

	if (!inherits(stp, "stp"))
		stop("'stp' must be an object of class 'stp'")

	#RGAS <- 8.314510       # J mol-1 deg-1 (perfect Gas)
	R <- 83.14472           # mol bar deg-1
	TK <- .TK(stp$T)		# T in degrees Kelvin

	# Note: there is an error in Table 9 of Millero, 1995.
	# The coefficients -b0 and b1
	# have to be multiplied by 1.e-3!

	# There are some more errors!
	# the signs (+,-) of coefficients in Millero 95 do not
	# agree with Millero 79
	a0 <- c(-25.5, -15.82, -29.48, -25.60, -18.03, -9.78, -48.76, -46,
		-14.51, -23.12, -26.57, -29.48, -26.43, -14.80)
	a1 <- c(0.1271, -0.0219, 0.1622, 0.2324, 0.0466, -0.0090, 0.5304,
		0.5304, 0.1211, 0.1758, 0.2020, 0.1622, 0.0889, 0.002)
	a2 <- c(0.0, 0.0, 2.608*1e-3, -3.6246*1e-3, 0.316*1e-3, -0.942*1e-3,
		0.0, 0.0, -0.321*1e-3, -2.647*1e-3, -3.042*1e-3, 2.6080e-3, -0.000905,
		-0.4e-3)
	b0 <- c(-3.08*1e-3, 1.13*1e-3, -2.84*1e-3, -5.13*1e-3, -4.53*1e-3,
		-3.91*1e-3, -11.76*1e-3, -11.76*1e-3, -2.67*1e-3, -5.15*1e-3,
		-4.08*1e-3, -2.84e-3, -5.03e-3, 2.89e-3)
	b1 <- c(0.0877*1e-3, -0.1475*1e-3, 0.0, 0.0794*1e-3, 0.09*1e-3,
		0.054*1e-3, 0.3692*1e-3, 0.3692*1e-3, 0.0427*1e-3, 0.09*1e-3,
		0.0714*1e-3, 0.0, 0.0814e-3, 0.054e-3)
	b2 <- c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

	deltav <- a0[i] + a1[i] * stp$T + a2[i] * stp$T^2
	deltak <- b0[i] + b1[i] * stp$T + b2[i] * stp$T^2
	lnkpok0 <- -(deltav / (R * TK)) * stp$P + (0.5 * deltak / (R * TK)) * stp$P^2
	return(exp(lnkpok0))
}

".checkUnit" <- function (unit)
{
	# Default value
	if (is.null(unit))
		unit <- "mol/kg-soln" else unit <- as.character(unit[1])

	# Separate prefix from the unit
	pre <- sub("mol.*$", "", unit)
	if (pre == unit) {
		pre <- substr(unit, 1, nchar(unit) - 1) # if 'M', or 'm'
		unit <- substr(unit, nchar(unit), 100)
	} else unit <- sub("^.*mol", "mol", unit)

	# Check prefix, can be nothing, m/milli, or mu/micro, n/nano
	if (pre %in% c("m", "milli")) {
		pre <- "m"
	} else if (pre %in% c("u", "micro")) {
		pre <- "u"
	} else if (pre %in% c("n", "nano")) {
		pre <- "n"
	} else pre <- ""

	# Check unit
	if (!unit %in% c("mol/kg-soln", "mol/kg-H2O", "mol/L", "molarity",
		"molality", "molinity", "molar", "molal", "molin", "M", "m", "i"))
		stop("'unit' must be of the valid seacarb units")

	# Convert from synonyms
	if (unit %in% c("molarity", "molar", "M")) unit <- "mol/L"
	if (unit %in% c("molality", "molal", "m")) unit <- "mol/kg-H2O"
	if (unit %in% c("molinity", "molin", "i")) unit <- "mol/kg-soln"

	return(paste(pre, unit, sep = ""))
}

".checkpHscale" <- function (pHscale)
{
	# Set default or force coercion of pHscale
	if (is.null(pHscale))
		pHscale <- "total scale" else pHscale <- as.character(pHscale[1])

	# Check value and convert value
	pHscale <- switch(pHscale,
		"free scale" = "free scale",
		"free" = "free scale",
		"f" = "free scale",
		"F" = "free scale",
		"total scale" = "total scale",
		"total" = "total scale",
		"t" = "total scale",
		"T" = "total scale",
		"seawater scale" = "seawater scale",
		"seawater" = "seawater scale",
		"sws" = "seawater scale",
		"SWS" = "seawater scale",
		stop("'pHscale' must be 'free scale', 'total scale' or 'seawater scale'")
	)
	return(pHscale)
}

".convUnit" <- function (fromUnit, toUnit, stp)
{
	# Check units
	fromUnit <- .checkUnit(fromUnit)
	toUnit <- .checkUnit(toUnit)
	# and stp
	if (is.null(stp) || !inherits(stp, "stp"))
		stop("Conversion requires correct definition of 'stp' conditions")

	# Calculation of conversion factor from one unit to the other, given STP
	conv <- paste(fromUnit, "->", toUnit)

	# Pre(fix) is 'm' for milli-, 'u' for micro-, 'n' for nano-
	# then it is stripped out
	fromPre <- sub("mol/.*$", "", fromUnit)
	fromUnit <- sub("^(mmol|umol|nmol)", "mol", fromUnit)
	toPre <- sub("mol/.*$", "", toUnit)
	toUnit <- sub("^(mmol|umol|nmol)", "mol", toUnit)

	# The calculated conversion factor: [toUnit] = res * [fromUnit]
	res <- 1

	# 1) If there is 'mmol/X', 'umol/X' or 'nmol/X', transform into 'mol/X'
	if (fromPre == "m") {
		res <- res / 1000
	} else if (fromPre == "u") {
		res <- res / 1000000
	} else if (fromPre == "n") {
		res <- res / 1000000000
	}

	# 2) If fromUnit is molality, convert to molinity
	# We use the following conversion formula from Dickson & Riley (1979):
	# molinity (mol/kg-soln) = (1 - 0.001005 * S) * molarity (mol/kg-H2O)
	if (fromUnit == "mol/kg-H2O") {
		res <- res * (1 - 0.001005 * stp$S)
		fromUnit <- "mol/kg-soln"
	}

	# 3) If toUnit is molarity, transform molinity into molarity, using rho
	if (toUnit == "mol/L") res <- res * rho(stp) / 1000

	# 4) If fromUnit is molarity, transform to molinity using rho
	if (fromUnit == "mol/L") {
		res <- res / rho(stp) * 1000 	# rho is in kg/m^3 => divide it by 1000
		fromUnit <- "mol/kg-soln"
	}

	# 5) If toUnit is molality, use inverse transformation (as 2) from molinity
	if (toUnit == "mol/kg-H2O") res <- res / (1 - 0.001005 * stp$S)

	# 6) Apply toPre ('m' for milli-, 'u' for micro-, 'n' for nano-)
	if (toPre == "m") {
		res <- res * 1000
	} else if (toPre == "u") {
		res <- res * 1000000
	} else if (toPre == "n") {
		res <- res * 1000000000
	}

	attr(res, "conv") <- conv
	attr(res, "stp") <- stp
	return(res)
}

".compareStp" <- function (stp1, stp2, tolerance = c(S = 0.1, T = 0.1, P = 0.1))
{
	# Special case: if stp2 == NULL, always return TRUE!
	# because we consider using only stp1
	if (is.null(stp2)) return(TRUE)

	# Check arguments (we assume stp1 and stp2 are correct - internal use only)
	if (!is.numeric(tolerance) || length(tolerance) != 3 || any(tolerance < 0))
		stop("'tolerance' must be a vector of 3 positives numbers")

	Sequal <- isTRUE(all.equal(stp1$S, stp2$S,
		scale = 1, tolerance = tolerance[1]))
	Tequal <- isTRUE(all.equal(stp1$T, stp2$T,
		scale = 1, tolerance = tolerance[2]))
	Pequal <- isTRUE(all.equal(stp1$P, stp2$P,
		scale = 1, tolerance = tolerance[3]))

	return(Sequal & Tequal & Pequal)
}

".shortUnit" <- function (unit)
{
	res <- switch(unit,
		"mol/L" = "M",
		"mmol/L" = "mM",
		"umol/L" = "uM",
		"nmol/L" = "nM",
		"mol/kg-H2O" = "m",
		"mmol/kg-H2O" = "mm",
		"umol/kg-H2O" = "um",
		"nmol/kg-H2O" = "nm",
		"mol/kg-soln" = "i",
		"mmol/kg-soln" = "mi",
		"umol/kg-soln" = "ui",
		"nmol/kg-soln" = "ni",
		stop("Unknown unit ", unit)
	)
	return(paste("_", res, sep = ""))
}

".specNames" <- function (label, unit = NULL)
{
	# Get short unit ID
	if (is.null(unit)) sunit <- "" else
		sunit <- .shortUnit(.checkUnit(unit))

	# 'label' is not defined? Return generic names
	if (is.null(label)) return(paste("C", 1:4, sunit, sep = ""))

	# Get species according to substance name given in 'label'
	res <- switch(as.character(label[1]),
		CT = c("CO2", "HCO3", "CO3", "C4"),
		BT = c("BOH3", "BOH4", "C3", "C4"),
		ST = c("HSO4", "SO4", "C3", "C4"),
		HST = c("H2S", "HS", "C3", "C4"),
		FT =  c("HF", "F", "C3", "C4"),
		SiT = c("SiOH4", "SiOOH3", "C3", "C4"),
		PT = c("H3PO4", "H2PO4", "HPO4", "PO4"),
		NHT = c("NH4", "NH3", "C3", "C4"), # Because NT also includes NO2/NO3!
		c("C1", "C2", "C3", "C4")
	)

	# Add short unit
	return(paste(res, sunit, sep = ""))
}

".prettyUnit" <- function (unit)
{
	return(switch(unit,
		"mol/L" = quote(mol.L^-1),
		"mmol/L" = quote(mmol.L^-1),
		"umol/L" = quote(mu*mol.L^-1),
		"nmol/L" = quote(nmol.L^-1),
		"mol/kg-soln" = quote(mol.kg-soln^-1),
		"mmol/kg-soln" = quote(mmol.kg-soln^-1),
		"umol/kg-soln" = quote(mu*mol.kg-soln^-1),
		"nmol/kg-soln" = quote(nmol.kg-soln^-1),
		"mol/kg-H2O" = quote(mol.kg-H[2]*O^-1),
		"mmol/kg-H2O" = quote(mmol.kg-H[2]*O^-1),
		"umol/kg-H2O" = quote(mu*mol.kg-H[2]*O^-1),
		"umol/kg-H2O" = quote(nmol.kg-H[2]*O^-1),
		unit
	))
}

".prettyLabel" <- function (label)
{
	# Strip unit and trailing identifiers out
	label <- sub("_.*$", "", label)

	# The table of conversion
	conv <- c(
		rho = quote(rho),
		H2O = quote(H[2]*O),
		H = quote(H^"+"),
		pH = quote(pH),
		pHT = quote(pH[T]),
		pHF = quote(pH[F]),
		pHSWS = quote(pH[SWS]),
		pHNBS = quote(pH~~(NBS)),
		pHNIST = quote(pH~~(NIST)),
		OH = quote(OH^"-"),
		SW = quote(seawater),
		Cl = quote(Cl^"-"),
		Br = quote(Br^"-"),
		Na = quote(Na^"+"),
		K = quote(K^"+"),
		Ca = quote(Ca^"2+"),
		Mg = quote(Mg^"2+"),
		Sr = quote(Sr^"2+"),
		CaCO3 = quote(CaCO[3]),
		O2 = quote(O[2]),
		N2 = quote(N[2]),
		BT = quote(B[T]),
		BOH4 = quote(B(OH)[4]^"-"),
		BOH3 = quote(B(OH)[3]),
		CT = quote(C[T]),
		AT = quote(A[T]),
		CO2 = quote(CO[2]^"*"),
		fCO2 = quote(italic(f)(CO[2])),
		HCO3 = quote(HCO[3]^"-"),
		CO3 = quote(CO[3]^"2-"),
		FT = quote(F[T]),
		F = quote(F^"-"),
		HF = quote(HF),
		NT = quote(N[T]),
		NHT = quote(NH[3] + NH[4]^"+"),
		NH3 = quote(NH[3]),
		NH4 = quote(NH[4]^"+"),
		NO2 = quote(NO[2]^"-"),
		NO3 = quote(NO[3]^"-"),
		PT = quote(P[T]),
		PO4 = quote(PO[4]^"3-"),
		HPO4 = quote(HPO[4]^"2-"),
		H2PO4 = quote(H[2]*PO[4]^"-"),
		H3PO4 = quote(H[3]*PO[4]),
		SiT = quote(Si[T]),
		SiOOH3 = quote(SiO(OH)[3]^"-"),
		SiOH4 = quote(Si(OH)[4]),
		ST = quote(S[T]),
		SO4 = quote(SO[4]^"2-"),
		HSO4 = quote(HSO[4]^"-")
	)
	# Make sure that all labels are present in the table
	missing <- !(label %in% names(conv))
	if(any(missing)) {
		conv2 <- label[missing]
		names(conv2) <- conv2
		conv <- c(conv, conv2)
	}
	# Do the substitution
	res <- conv[label]
	return(as.expression(res))
}
