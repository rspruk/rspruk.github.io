* test_easc.do -- verifies the easc command against the paper's Python results
* Run from the folder containing easc.ado, easc.sthlp and iran_easc.csv:
*     do test_easc.do
clear all
set more off
adopath + "`c(pwd)'"

import delimited iran_easc.csv, clear
label var lgdppc  "log GDP per capita, PPP"
label var oilrent "oil rents, % of GDP, 1996-2011 mean"
label var lbrent  "log Brent price"

* ---------- 1. Exposure-feasible pool (21 donors): the paper's main specification
easc lgdppc, id(id) time(year) treated(1) trperiod(2012) exposure(oilrent) shock(lbrent) ///
    phases(2012 2015 2016 2017 2018 2024) placebo graph(all) graphprefix(iran_)

* Expected (Table 5 and Table 4 of the paper, exposure-feasible pool):
*   EI(w_SC) = 15.527   d_min = 0.000   RMSPE SC = 0.0144   RMSPE EBSC-min = 0.0526   R2 = 0.915
*   phase      SC        EASC      EBSC-min   100*phi   kappa   kappaEA
*   2012-2015  -0.1271   -0.1222   -0.1089    -0.031    12.12   12.11
*   2016-2017  -0.0873   -0.0107   -0.0218    -0.493    28.12   27.54
*   2018-2024  -0.0894    0.1094    0.0536    -1.280    40.40   39.95
*   SC   2012-15: p = 0.045, noise sd = 0.048, CI(M=0) = [-0.221, -0.033], CI(M=0.5) = [-0.271, 0.017], M* = 0.27
*   EASC 2012-15: p = 0.091, noise sd = 0.052, CI(M=0) = [-0.224, -0.021], M* = 0.21
*   weights: Mexico 0.369, Nigeria 0.254, Thailand 0.243, Kazakhstan 0.133 (ids 5, 3, 13, 2)
matrix list e(phases)
matrix list e(w)
matrix list e(diag)

* ---------- 2. Normaliser pool only (19 donors): Iran outside the exposure hull
easc lgdppc if normaliser == 1, id(id) time(year) treated(1) trperiod(2012) exposure(oilrent) shock(lbrent) ///
    phases(2012 2015 2016 2017 2018 2024)

* Expected: d_min = 6.313, EI = 15.527, RMSPE EBSC-min = 0.3486 (balance infeasible), R2 = 0.921
*   2012-2015  SC -0.1271  EASC -0.1515  EBSC-min -0.7007   kappa 7.86  kappaEA 7.67
*   2016-2017  SC -0.0873  EASC -0.0328  EBSC-min -0.6631   kappa 18.95 kappaEA 17.31
*   2018-2024  SC -0.0894  EASC  0.0862  EBSC-min -0.7375   kappa 26.72 kappaEA 25.29

* ---------- 3. Save the annual gaps and plot them
easc lgdppc, id(id) time(year) treated(1) trperiod(2012) exposure(oilrent) ///
    phases(2012 2015 2016 2017 2018 2024) saving(iran_gaps.dta) replace
preserve
use iran_gaps.dta, clear
twoway (line gap_sc year, lpattern(dash)) (line gap_easc year), yline(0) xline(2011.5) ///
    legend(order(1 "Synthetic control" 2 "Exposure-augmented SC")) ytitle("Gap, log points")
restore
