import delimited simulation.tab, case(preserve) 
tsset T
gen valid = 0
replace valid = 1 if d == 0 & L.d == 0 & L2.d == 0 & L3.d == 0
gen durPr = 1 + 9 * (bL / (bS + bL))
gen nx = gdp - c
gen nxy = nx / gdp
gen spS = 1.032 * (1 / qS - 1)
gen spL = (0.0712 + 0.032) * (1 / qL - 1)
gen ellL = bL - (1-0.0712) * L.bL
replace ellL = 0.0 if ellL < 0.0
gen durIssue = 1 + 9 * (ellL / (bS + ellL))
compress
