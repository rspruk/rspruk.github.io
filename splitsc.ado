*! splitsc 0.9.6  Time-split honest inference for synthetic controls with imperfect fit
*! Implements the trend-balancing time-split interval, the curvature class with the slope
*! unrestricted, the breakdown curvature, placebo curvatures, the breakdown persistence under
*! AR(1) dependence and the placebo-calibrated drift bound (Theorem 3) of the accompanying paper.
program define splitsc, rclass
    version 15
    syntax varname(numeric) [if] [in], Unit(varname) Time(varname numeric) TReated(string) ///
        TRPeriod(real) TRAINend(real) [ NOINTercept EXClude(numlist) Curvature(real 0) ///
        Level(cilevel) PHIbar(numlist >=0 <1) REPs(integer 20000) SEED(string) EXCH NOPLacebo ///
        STep(real 0.25) GRaph(string) GName(string) DEBUG ]

    marksample touse, novarlist
    markout `touse' `unit' `time', strok
    if "`seed'" != "" set seed `seed'
    if `trainend' >= `trperiod' {
        di as err "trainend() must be earlier than trperiod()"
        exit 198
    }
    if `step' <= 0 | `reps' < 100 {
        di as err "step() must be positive and reps() at least 100"
        exit 198
    }
    if `curvature' < 0 {
        di as err "curvature() must be nonnegative"
        exit 198
    }
    if "`phibar'" == "" local phibar 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9

    * identify the treated unit
    tempvar uid
    qui egen long `uid' = group(`unit') if `touse'
    capture confirm string variable `unit'
    if !_rc qui su `uid' if `unit' == `"`treated'"' & `touse', meanonly
    else {
        capture confirm number `treated'
        if _rc {
            di as err "treated() must be a number when `unit' is numeric"
            exit 198
        }
        qui su `uid' if float(`unit') == float(`treated') & `touse', meanonly
    }
    if r(N) == 0 {
        di as err `"treated unit `treated' not found in `unit'"'
        exit 198
    }
    local trid = r(min)
    qui levelsof `unit' if `touse' & `uid' != `trid', local(donors)

    * the exclusion list is passed to Mata as a row vector
    tempname EXM PHM
    if "`exclude'" != "" {
        local nx : word count `exclude'
        matrix `EXM' = J(1, `nx', .)
        forvalues k = 1/`nx' {
            matrix `EXM'[1, `k'] = `: word `k' of `exclude''
        }
    }
    else matrix `EXM' = J(1, 1, .)
    local np : word count `phibar'
    matrix `PHM' = J(1, `np', .)
    forvalues k = 1/`np' {
        matrix `PHM'[1, `k'] = `: word `k' of `phibar''
    }

    local icpt = ("`nointercept'" == "")
    local doexch = ("`exch'" != "")
    local doplac = ("`noplacebo'" == "")
    local dbg = ("`debug'" != "")
    if `doexch' & !`doplac' {
        di as err "option exch requires the placebo pool; do not combine it with noplacebo"
        exit 198
    }

    * figures
    local graph = lower(strtrim(`"`graph'"'))
    if "`graph'" == "all" {
        local graph "paths gap sensitivity persistence"
        if `doplac' local graph "`graph' placebo"
    }
    local glist
    foreach g of local graph {
        if inlist("`g'", "path", "paths") local g paths
        else if "`g'" == "gap" local g gap
        else if inlist("`g'", "sens", "sensitivity", "curvature") local g sensitivity
        else if inlist("`g'", "plac", "placebo", "placebos") local g placebo
        else if inlist("`g'", "pers", "persistence", "phibar") local g persistence
        else {
            di as err "graph(): `g' is not a figure. Use paths, gap, sensitivity, placebo, persistence or all"
            exit 198
        }
        if "`g'" == "placebo" & !`doplac' {
            di as err "graph(placebo) needs the placebo pool; do not combine it with noplacebo"
            exit 198
        }
        local glist `glist' `g'
    }
    local glist : list uniq glist
    local dofull : list posof "persistence" in glist
    local dofull = (`dofull' > 0)
    if "`gname'" == "" local gname splitsc_
    if "`glist'" != "" {
        capture confirm name `gname'x
        if _rc {
            di as err "gname() must be a valid name prefix"
            exit 198
        }
    }

    mata: _splitsc_main("`varlist'", "`uid'", "`time'", "`touse'", `trid', `trperiod', `trainend', ///
        "`EXM'", `icpt', `curvature', `level', "`PHM'", `reps', `doexch', `doplac', `step', `dbg', `dofull')

    * ---------------------------------------------------------------- display
    local lv = `level'
    di
    di as txt "Time-split synthetic control" _col(50) "Number of donors  = " as res %8.0f r(J)
    di as txt "Weights: " as res cond(`icpt', "simplex with intercept (de-meaned)", "simplex without intercept") ///
        _col(50) as txt "Training periods  = " as res %8.0f r(n_train)
    di as txt "Anchor: trend-balancing" _col(50) "Validation periods= " as res %8.0f r(n)
    di as txt _col(50) "Post periods      = " as res %8.0f r(m)
    di as txt _col(50) "Residual df       = " as res %8.0f r(nu)
    di as txt "{hline 78}"
    di as txt "Average effect" _col(22) "Estimate" _col(36) "Std. err." _col(52) "[`lv'% conf. interval]"
    di as txt "{hline 78}"
    di as txt "Linear class C=0" _col(18) as res %12.4g r(theta) _col(32) %12.4g r(se) _col(48) %12.4g r(lb) _col(62) %12.4g r(ub)
    if `curvature' > 0 di as txt "Curvature C=" as res %-8.4g `curvature' _col(32) as txt "" _col(48) as res %12.4g r(lb_C) _col(62) %12.4g r(ub_C)
    if `doplac' di as txt "Curvature C=C_don" _col(48) as res %12.4g r(lb_Cdon) _col(62) %12.4g r(ub_Cdon)
    if `doexch' di as txt "Theorem 3" _col(18) as res %12.4g r(exch_theta) _col(48) %12.4g r(exch_lb) _col(62) %12.4g r(exch_ub)
    di as txt "{hline 78}"
    di as txt "Curvature exposure kappa_C" _col(40) "= " as res %12.5g r(kappa)
    di as txt "Breakdown curvature C*" _col(40) "= " as res %12.5g r(Cstar)
    if `doplac' {
        di as txt "Largest placebo curvature C_don" _col(40) "= " as res %12.5g r(Cdon)
        di as txt "Placebos with curvature above C*" _col(40) "= " as res %12.0f r(above) as txt " of " as res r(J)
    }
    if r(phistar) < . di as txt "Breakdown persistence phibar*" _col(40) "= " as res %12.2f r(phistar)
    else             di as txt "Breakdown persistence phibar*" _col(40) "> " as res %12.2f r(phimax)
    if `doexch' {
        di as txt "Theorem 3: k_sim = " as res %6.3f r(ksim) as txt ", M_hat = " as res %10.4g r(Mhat) ///
            as txt ", h = " as res %8.4g r(exch_h) as txt ", guarantee >= " as res %5.3f r(exch_guarantee)
    }
    di as txt "{hline 78}"
    di as txt "Notes: the interval is valid in finite samples under Gaussian, serially independent"
    di as txt "noise and the stated bias class. Rows at C_don and the counts of placebos are a"
    di as txt "sensitivity analysis, since C_don uses donor outcomes after the training block."

    * donor weights above 0.001
    tempname W
    matrix `W' = r(weights)
    di as txt "Donor weights (above 0.001)"
    local k 0
    foreach d of local donors {
        local ++k
        if `W'[1, `k'] > 0.001 di as txt "  " `"`d'"' _col(40) as res %8.4f `W'[1, `k']
    }

    * ---------------------------------------------------------------- figures
    if "`glist'" != "" {
        tempname GM PM TM PLM PSM HOLD
        matrix `GM' = r(gap)
        matrix `PM' = r(paths)
        matrix `TM' = r(trend)
        if `doplac' matrix `PLM' = r(placebo)
        else matrix `PLM' = J(1, 2, .)
        if `dofull' matrix `PSM' = r(persist)
        else matrix `PSM' = J(1, 4, .)
        local cdon = cond(`doplac', r(Cdon), -1)
        local ylab : variable label `varlist'
        if `"`ylab'"' == "" local ylab `varlist'
        local tlab : variable label `time'
        if `"`tlab'"' == "" local tlab `time'
        local tfmt : format `time'
        local gargs theta(`=r(theta)') se(`=r(se)') tq(`=r(tcrit)') kappa(`=r(kappa)') ///
            cstar(`=r(Cstar)') cdon(`cdon') phistar(`=r(phistar)') cuser(`curvature') ///
            level(`level') trperiod(`trperiod')
        _return hold `HOLD'
        capture noisily _splitsc_graph, graphs(`glist') prefix(`gname') gm(`GM') pm(`PM') tm(`TM') ///
            plm(`PLM') psm(`PSM') `gargs' treated(`"`treated'"') donors(`"`donors'"') ///
            ylab(`"`ylab'"') tlab(`"`tlab'"') tfmt(`tfmt')
        local grc = _rc
        _return restore `HOLD'
        if `grc' exit `grc'
    }

    return add
    return local donors `"`donors'"'
    return local treated `"`treated'"'
    return local cmd "splitsc"
end

* ---------------------------------------------------------------- figures
program define _splitsc_graph
    version 15
    syntax , graphs(string) prefix(string) gm(name) pm(name) tm(name) plm(name) psm(name) ///
        theta(real) se(real) tq(real) kappa(real) cstar(real) cdon(real) phistar(string) ///
        cuser(real) level(real) trperiod(real) treated(string) donors(string) ///
        ylab(string) tlab(string) tfmt(string)

    local hw = `tq' * `se'
    local lvl = strofreal(`level')
    preserve
    quietly {
        drop _all
        local nobs = max(rowsof(`gm'), 201, rowsof(`plm'), rowsof(`psm'))
        set obs `nobs'
        svmat double `gm', names(col)
        svmat double `pm', names(col)
        format time `tfmt'
        gen long _seq = _n
        gen byte _blk = cond(train == 1, 1, cond(validation == 1, 2, cond(post == 1, 3, 0))) if time < .

        * period spacing and block edges
        tempvar dt
        gen double `dt' = time[_n + 1] - time if time < . & time[_n + 1] < .
        su `dt', meanonly
        local d2 = cond(r(min) < ., r(min) / 2, 0.5)
        su time if _blk == 2, meanonly
        local va = r(min) - `d2'
        local vb = r(max) + `d2'
        local x0 = r(max)
        su time if _blk == 3, meanonly
        local pa = r(min) - `d2'
        local pb = r(max) + `d2'
        local xbar = r(mean)
        su time if _blk == 1, meanonly
        local ta = r(min)
        local tb = r(max)
        su gap if _blk == 3, meanonly
        local mq = r(mean)
        local b0 = `tm'[1, 1]
        local b1 = `tm'[1, 2]
        local tbar = `b0' + `b1' * `xbar'
    }
    local shade "lwidth(none)"
    local xl `"xtitle(`"`tlab'"') xlabel(, format(`tfmt'))"'
    local common `"graphregion(color(white)) plotregion(color(white)) legend(region(lcolor(white)) size(small) position(6) symxsize(*0.6) keygap(*0.8) colgap(*1.5))"'

    foreach g of local graphs {
        * ---------------------------------------------------------- paths
        if "`g'" == "paths" {
            quietly su treated if time < ., meanonly
            local lo = r(min)
            local hi = r(max)
            quietly su synthetic if time < ., meanonly
            local lo = min(`lo', r(min))
            local hi = max(`hi', r(max))
            local pad = (`hi' - `lo') * 0.05
            local lo = `lo' - `pad'
            local hi = `hi' + `pad'
            twoway (function y = `hi', range(`va' `vb') recast(area) base(`lo') color(navy%12) `shade') ///
                (function y = `hi', range(`pa' `pb') recast(area) base(`lo') color(cranberry%10) `shade') ///
                (line synthetic time, lcolor(navy) lpattern(dash) lwidth(medthick)) ///
                (line treated time, lcolor(black) lwidth(medthick)), ///
                yscale(range(`lo' `hi')) ytitle(`"`ylab'"') `xl' ///
                legend(order(4 `"`treated'"' 3 "Synthetic control" 1 "Validation" 2 "Post-treatment") rows(1)) ///
                title("Treated unit and synthetic control", size(medium)) ///
                note("Weights fitted on the training block only (unshaded, before the validation block).", size(vsmall)) ///
                `common' name(`prefix'paths, replace)
        }
        * ---------------------------------------------------------- gap
        if "`g'" == "gap" {
            quietly {
                capture drop _cx _clo _chi _cmid _dlo _dhi
                gen double _cx = `xbar' in 1
                gen double _clo = `tbar' - `hw' in 1
                gen double _chi = `tbar' + `hw' in 1
                gen double _cmid = `mq' in 1
                local extra
                local dleg
                if `cdon' >= 0 {
                    gen double _dlo = _clo - `cdon' * `kappa' in 1
                    gen double _dhi = _chi + `cdon' * `kappa' in 1
                    local extra (rcap _dlo _dhi _cx, lcolor(gs9) lwidth(medium) msize(large))
                    local dleg 8 "Band at C = C_don"
                }
                su gap if time < ., meanonly
                local lo = min(r(min), _clo[1], 0)
                local hi = max(r(max), _chi[1], 0)
                if `cdon' >= 0 {
                    local lo = min(`lo', _dlo[1])
                    local hi = max(`hi', _dhi[1])
                }
                local f0 = min(`b0' + `b1' * (`va' + `d2'), `b0' + `b1' * (`pb' - `d2'))
                local f1 = max(`b0' + `b1' * (`va' + `d2'), `b0' + `b1' * (`pb' - `d2'))
                local lo = min(`lo', `f0')
                local hi = max(`hi', `f1')
                local pad = (`hi' - `lo') * 0.05
                local lo = `lo' - `pad'
                local hi = `hi' + `pad'
                local fa = `va' + `d2'
                local fb = `pb' - `d2'
                local qa = `pa' + `d2'
            }
            twoway (function y = `hi', range(`va' `vb') recast(area) base(`lo') color(navy%12) `shade') ///
                (function y = `hi', range(`pa' `pb') recast(area) base(`lo') color(cranberry%10) `shade') ///
                (line gap time, lcolor(black) lwidth(medium)) ///
                (scatter gap time if _blk > 0, mcolor(black) msymbol(O) msize(small)) ///
                (scatter gap time if _blk == 0, mcolor(black) msymbol(Oh) msize(small)) ///
                (function y = `b0' + `b1' * x, range(`fa' `fb') lcolor(navy) lpattern(dash) lwidth(medthick)) ///
                (function y = `mq', range(`qa' `fb') lcolor(cranberry) lwidth(medthick)) ///
                `extra' ///
                (rcap _clo _chi _cx, lcolor(cranberry) lwidth(thick) msize(large)) ///
                (scatter _cmid _cx, mcolor(cranberry) msymbol(D) msize(medlarge)), ///
                yline(0, lcolor(gs10)) yscale(range(`lo' `hi')) ytitle("Gap (treated minus synthetic)") `xl' ///
                legend(order(6 "Validation trend, extrapolated" 7 "Post-treatment mean" ///
                    `=cond(`cdon' >= 0, 9, 8)' "`lvl'% band for the counterfactual, linear class" `dleg') cols(2)) ///
                title("Gap and the trend-balancing estimate", size(medium)) ///
                note("Diamond: post-treatment mean gap. Bands are centred on the validation trend at the mean post-treatment period." ///
                     "The effect estimate is the distance from the diamond to the centre of the band. The interval excludes zero when the" ///
                     "diamond lies outside the band. Shaded: validation and post-treatment blocks. Hollow markers: excluded periods.", size(vsmall)) ///
                `common' name(`prefix'gap, replace)
        }
        * ---------------------------------------------------------- sensitivity to curvature
        if "`g'" == "sensitivity" {
            quietly {
                capture drop _cg _slo _shi _sth _pc _py
                local cmax = max(1.25 * max(`cdon', 0), 2 * `cstar', 1.25 * `cuser')
                if `cmax' <= 0 local cmax = 1
                if _N < 201 set obs 201
                gen double _cg = `cmax' * (_n - 1) / 200 in 1/201
                gen double _slo = `theta' - `hw' - _cg * `kappa'
                gen double _shi = `theta' + `hw' + _cg * `kappa'
                gen double _sth = `theta' if _cg < .
                su _slo, meanonly
                local lo = min(r(min), 0)
                su _shi, meanonly
                local hi = max(r(max), 0)
                local pad = (`hi' - `lo') * 0.06
                local lo = `lo' - `pad'
                local hi = `hi' + `pad'
                local ptick
                local pleg
                local xls
                if `cdon' >= 0 {
                    if _N < rowsof(`plm') set obs `=rowsof(`plm')'
                gen double _pc = `plm'[_n, 2] in 1/`=rowsof(`plm')'
                    gen double _py = `lo' if _pc < .
                    local ptick (scatter _py _pc, msymbol(pipe) msize(huge) mcolor(gs6))
                    local pleg 4 "Placebo curvatures"
                    local xls xline(`cdon', lcolor(gs8) lpattern(dash))
                }
                local tpos = `hi' - `pad'
                local tdon
                if `cuser' > 0 local xls `xls' xline(`cuser', lcolor(navy) lpattern(shortdash))
                if `cuser' > 0 local tdon text(`tpos' `cuser' " C", place(e) color(navy) size(small))
                if `cdon' >= 0 local tdon `tdon' text(`tpos' `cdon' " C_don", place(e) color(gs6) size(small))
            }
            twoway (rarea _slo _shi _cg, color(navy%25) lwidth(none)) ///
                (line _sth _cg, lcolor(navy) lwidth(medthick)) ///
                (function y = 0, range(0 `cmax') lcolor(gs10)) ///
                `ptick', ///
                xline(`cstar', lcolor(cranberry) lwidth(medthick)) `xls' ///
                text(`tpos' `cstar' " C*", place(e) color(cranberry) size(small)) ///
                `tdon' ///
                yscale(range(`lo' `hi')) ytitle("Average effect") xtitle("Curvature bound C") ///
                legend(order(1 "`lvl'% interval" 2 "Estimate" `pleg') rows(1)) ///
                title("Sensitivity to the curvature bound", size(medium)) ///
                note("C*: breakdown curvature, where the interval first includes zero. C_don: largest placebo curvature. C: the bound in curvature().", size(vsmall)) ///
                `common' name(`prefix'sensitivity, replace)
        }
        * ---------------------------------------------------------- placebo curvatures
        if "`g'" == "placebo" {
            quietly {
                capture drop _pcv _don _rk
                capture label drop _ssc_dl
                local J = rowsof(`plm')
                if _N < `J' set obs `J'
                gen double _pcv = `plm'[_n, 2] in 1/`J'
                gen strL _don = ""
                local k 0
                foreach d of local donors {
                    local ++k
                    replace _don = `"`d'"' in `k'
                }
                gsort -_pcv
                gen int _rk = _n if _pcv < .
                forvalues i = 1/`J' {
                    label define _ssc_dl `i' `"`=_don[`i']'"', add
                }
                label values _rk _ssc_dl
                count if _pcv > `cstar' & _pcv < .
                local nab = r(N)
                local ys = min(max(4.5, 0.22 * `J'), 9)
                local ls = cond(`J' > 25, "tiny", cond(`J' > 15, "vsmall", "small"))
            }
            twoway (scatter _rk _pcv if _pcv <= `cstar', mcolor(navy) msymbol(O)) ///
                (scatter _rk _pcv if _pcv > `cstar' & _pcv < ., mcolor(cranberry) msymbol(D)), ///
                xline(`cstar', lcolor(cranberry) lwidth(medthick)) ///
                ylabel(1/`J', valuelabel angle(0) labsize(`ls') nogrid) yscale(reverse) ytitle("") ///
                xtitle("Placebo curvature") ///
                legend(order(1 "Below C*" 2 "Above C*") rows(1)) ///
                title("Placebo curvatures and the breakdown curvature", size(medium)) ///
                note("`nab' of `J' placebos have curvature above C* (vertical line).", size(vsmall)) ///
                ysize(`ys') `common' name(`prefix'placebo, replace)
            quietly sort _seq
        }
        * ---------------------------------------------------------- persistence
        if "`g'" == "persistence" {
            quietly {
                capture drop _ph _plo _phi _pth
                local np = rowsof(`psm')
                if _N < `np' set obs `np'
                gen double _ph = `psm'[_n, 1] in 1/`np'
                gen double _plo = `psm'[_n, 3] in 1/`np'
                gen double _phi = `psm'[_n, 4] in 1/`np'
                gen double _pth = `theta' if _ph < .
                su _plo, meanonly
                local lo = min(r(min), 0)
                su _phi, meanonly
                local hi = max(r(max), 0)
                su _ph, meanonly
                local pmx = r(max)
                local xls
                local pnote "The interval excludes zero over the whole grid."
                if real("`phistar'") < . {
                    local xls xline(`phistar', lcolor(cranberry) lwidth(medthick))
                    local pnote "Vertical line: breakdown persistence phibar*, where the interval first includes zero."
                }
            }
            twoway (rarea _plo _phi _ph, color(navy%25) lwidth(none)) ///
                (rcap _plo _phi _ph, lcolor(navy) msize(medium)) ///
                (line _pth _ph, lcolor(navy) lwidth(medthick)) ///
                (function y = 0, range(0 `pmx') lcolor(gs10)), ///
                `xls' ytitle("Average effect") xtitle("Largest AR(1) coefficient phibar") ///
                legend(order(1 "`lvl'% interval" 3 "Estimate") rows(1)) ///
                title("Sensitivity to serial dependence", size(medium)) ///
                note("`pnote'", size(vsmall)) ///
                `common' name(`prefix'persistence, replace)
        }
    }
    restore
end

version 15
capture mata: mata drop _ssc_nnls() _ssc_simplex() _ssc_K() _ssc_kappa() _ssc_rank() _ssc_knum() _splitsc_main()
mata:
mata set matastrict off

// ---------------------------------------------------------------- Lawson-Hanson nonnegative least squares
real colvector _ssc_nnls(real matrix A, real colvector b)
{
    n = cols(A); x = J(n, 1, 0); P = J(n, 1, 0); z = J(n, 1, 0)
    tol = 1e-10 * max((1, max(abs(A))))
    w = A' * (b - A * x); it = 0
    while (sum(P :== 0) > 0 & it < 3 * n) {
        wz = w :* (P :== 0) + (-1e300) :* P
        if (max(wz) <= tol) break
        ii = .; ww = .
        maxindex(wz, 1, ii, ww)
        P[ii[1]] = 1
        for (inner = 1; inner <= 3 * n; inner++) {
            idx = selectindex(P)
            z = J(n, 1, 0)
            z[idx] = pinv(A[., idx]) * b
            if (min(z[idx]) > 0) break
            zi = z[idx]
            neg = select(idx, zi :<= 0)
            alpha = min(x[neg] :/ (x[neg] - z[neg]))
            if (alpha >= .) alpha = 0
            x = x + alpha * (z - x)
            P = P :* (x :> tol)
            x = x :* P
        }
        x = z
        w = A' * (b - A * x)
        it++
    }
    return(x)
}

// simplex weights: min ||b - cc - A w||, w >= 0, sum(w) = 1, cc free if intercept
real colvector _ssc_simplex(real matrix A, real colvector b, real scalar intercept, real scalar cc)
{
    if (intercept) {
        Ac = A :- mean(A); bc = b :- mean(b)
    }
    else {
        Ac = A; bc = b
    }
    s = max((max(abs(Ac)), max(abs(bc)), 1e-12)) * 1e4
    w = _ssc_nnls(Ac \ J(1, cols(A), s), bc \ s)
    if (sum(w) > 0) w = w / sum(w)
    else w = J(cols(A), 1, 1 / cols(A))
    if (intercept) cc = mean(b - A * w)
    else cc = 0
    return(w)
}

// ---------------------------------------------------------------- Peano kernel and curvature exposure
real matrix _ssc_K(real colvector x, real rowvector s, real scalar x0)
{
    G = J(rows(x), 1, 1) * (s :>= x0)
    D = x * J(1, cols(s), 1) - J(rows(x), 1, 1) * s
    return(G :* (D :* (D :> 0)) + (1 :- G) :* ((-D) :* (D :< 0)))
}

real scalar _ssc_kappa(real colvector lam, real colvector xP, real colvector xQ, real scalar step)
{
    lo = min((min(xP), min(xQ))); hi = max((max(xP), max(xQ))); x0 = max(xP)
    ns = floor((hi - lo) / step + 1e-9) + 1
    s = (lo :+ (0::(ns - 1)) * step)'
    KQ = mean(_ssc_K(xQ, s, x0))
    KP = _ssc_K(xP, s, x0)
    return(sum(abs(KQ - lam' * KP)) * step)
}

real scalar _ssc_rank(real matrix X)
{
    sv = svdsv(X)
    return(sum(sv :> 1e-8 * max(sv)))
}

// ---------------------------------------------------------------- AR(1) numerical critical value
real scalar _ssc_knum(real colvector lam, real colvector xP, real colvector xQ, real scalar phibar,
                      real scalar alpha, real scalar R)
{
    n = rows(xP); m = rows(xQ); X = (J(n, 1, 1), xP, lam)
    Hm = I(n) - X * pinv(X); nu = n - _ssc_rank(X); a = sqrt(1 / m + lam' * lam)
    tt = xP \ xQ; kmax = 0
    for (g = 0; g <= 10; g++) {
        phi = phibar * g / 10
        if (phi == 0) Rm = I(n + m)
        else Rm = phi :^ abs(tt * J(1, n + m, 1) - J(n + m, 1, 1) * tt')
        L = cholesky(Rm + 1e-12 * I(n + m))
        Z = rnormal(R, n + m, 0, 1) * L'
        xi = mean(Z[., (n + 1)..(n + m)]')' - Z[., 1..n] * lam
        E = Z[., 1..n] * Hm'
        sh = a * sqrt(rowsum(E :^ 2) / nu)
        rr = sort(abs(xi) :/ sh, 1)
        k = rr[ceil((1 - alpha) * R)]
        if (k > kmax) kmax = k
    }
    return(kmax)
}

// ---------------------------------------------------------------- main routine
void _splitsc_main(string scalar yv, string scalar uv, string scalar tv, string scalar touse,
                   real scalar trid, real scalar trperiod, real scalar trainend, string scalar exm,
                   real scalar intercept, real scalar C, real scalar level, string scalar phm,
                   real scalar R, real scalar doexch, real scalar doplac, real scalar step, real scalar dbg,
                   real scalar dofull)
{
    st_rclear()
    D = st_data(., (yv, uv, tv), touse)
    D = sort(D, (2, 3))
    ids = uniqrows(D[., 2]); times = uniqrows(D[., 3]); N = rows(ids); T = rows(times)
    if (rows(D) != N * T) {
        errprintf("the panel must be balanced: every unit needs one observation in every period\n")
        exit(459)
    }
    if (dbg) printf("[1] data %f x %f, units %f, periods %f\n", rows(D), cols(D), N, T)
    Ymat = colshape(D[., 1], T)
    if (colshape(D[., 3], T) != J(N, 1, times')) {
        errprintf("each unit must be observed once in every period\n")
        exit(459)
    }
    if (hasmissing(Ymat)) {
        errprintf("the outcome has missing values in the estimation sample\n")
        exit(416)
    }
    tr = selectindex(ids :== trid)
    order = tr \ select((1::N), ids :!= trid)
    Ymat = Ymat[order, .]
    nd = N - 1
    if (nd < 2) {
        errprintf("at least two donors are needed\n")
        exit(2001)
    }

    ex = st_matrix(exm)
    keep = J(T, 1, 1)
    for (k = 1; k <= cols(ex); k++) if (ex[k] < .) keep = keep :* (times :!= ex[k])
    pre = keep :* (times :< trperiod)
    trn = pre :* (times :<= trainend)
    val = pre :* (times :> trainend)
    pst = keep :* (times :>= trperiod)
    itr = selectindex(trn); iva = selectindex(val); ipo = selectindex(pst)
    n = rows(iva); m = rows(ipo); n1 = rows(itr)
    if (n < 3) {
        errprintf("at least three validation periods are needed (periods after trainend() and before trperiod())\n")
        exit(2001)
    }
    if (m < 1 | n1 < 2) {
        errprintf("not enough training or post-treatment periods\n")
        exit(2001)
    }
    xP = times[iva]; xQ = times[ipo]
    if (dbg) printf("[2] train %f, validation %f, post %f, xP %f x %f\n", n1, n, m, rows(xP), cols(xP))
    alpha = 1 - level / 100

    // weights on the training block and the gap
    y1 = Ymat[1, .]'
    Yc = Ymat[(2::N), .]
    if (dbg) printf("[2a] treated %f x %f, donors %f x %f\n", rows(y1), cols(y1), rows(Yc), cols(Yc))
    Atr = Yc[., itr]'
    btr = y1[itr]
    if (dbg) printf("[2b] training design %f x %f, target %f x %f\n", rows(Atr), cols(Atr), rows(btr), cols(btr))
    c0 = 0
    w = _ssc_simplex(Atr, btr, intercept, c0)
    if (dbg) printf("[2c] weights %f x %f, intercept %f x %f\n", rows(w), cols(w), rows(c0), cols(c0))
    fit = Yc' * w
    u = y1 - fit
    u = u :- c0
    uP = u[iva]; uQ = u[ipo]
    if (dbg) printf("[3] weights %f x %f, gap %f x %f, uP %f x %f, uQ %f x %f\n", rows(w), cols(w), rows(u), cols(u), rows(uP), cols(uP), rows(uQ), cols(uQ))

    // trend-balancing anchor and the Student interval
    A = (J(n, 1, 1), xP)
    lam = A * invsym(A' * A) * (1 \ mean(xQ))
    theta = mean(uQ) - lam' * uP
    X = (A, lam); nu = n - _ssc_rank(X)
    res = uP - X * pinv(X) * uP
    omega = sqrt(res' * res / nu); se = omega * sqrt(1 / m + lam' * lam)
    tq = invttail(nu, alpha / 2)
    if (dbg) printf("[4] anchor %f x %f, nu %f\n", rows(lam), cols(lam), nu)
    kap = _ssc_kappa(lam, xP, xQ, step)
    if (dbg) printf("[5] kappa %f\n", kap)
    Cstar = max((abs(theta) - tq * se, 0)) / kap

    st_numscalar("r(theta)", theta); st_numscalar("r(se)", se); st_numscalar("r(omega)", omega)
    st_numscalar("r(nu)", nu); st_numscalar("r(tcrit)", tq)
    st_numscalar("r(lb)", theta - tq * se); st_numscalar("r(ub)", theta + tq * se)
    st_numscalar("r(kappa)", kap); st_numscalar("r(Cstar)", Cstar)
    st_numscalar("r(lb_C)", theta - C * kap - tq * se); st_numscalar("r(ub_C)", theta + C * kap + tq * se)
    st_numscalar("r(C)", C); st_numscalar("r(level)", level)
    st_numscalar("r(J)", nd); st_numscalar("r(n)", n); st_numscalar("r(m)", m); st_numscalar("r(n_train)", n1)
    st_numscalar("r(intercept)", c0)
    st_matrix("r(weights)", w'); st_matrix("r(anchor)", lam')
    if (dbg) printf("[6] gap matrix parts %f %f %f %f %f rows\n", rows(times), rows(u), rows(trn), rows(val), rows(pst))
    st_matrix("r(gap)", (times, u, trn, val, pst))
    st_matrixcolstripe("r(gap)", (J(5, 1, ""), ("time" \ "gap" \ "train" \ "validation" \ "post")))
    st_matrix("r(paths)", (times, y1, fit :+ c0))
    st_matrixcolstripe("r(paths)", (J(3, 1, ""), ("ptime" \ "treated" \ "synthetic")))
    bt = invsym(A' * A) * A' * uP
    st_matrix("r(trend)", bt')
    st_matrixcolstripe("r(trend)", (J(2, 1, ""), ("intercept" \ "slope")))

    // placebos: each donor against the other donors, treated unit excluded
    Pi = J(0, 1, .)
    if (doplac) {
        xx = (xP \ xQ); xx = xx :- mean(xx); XQ = (J(rows(xx), 1, 1), xx, xx :^ 2)
        Pi = J(nd, 1, .); Cj = J(nd, 1, .)
        for (j = 1; j <= nd; j++) {
            oth = select((2::N), (2::N) :!= j + 1)
            yj = Ymat[j + 1, .]'
            Yo = Ymat[oth, .]
            cj0 = 0
            wj = _ssc_simplex(Yo[., itr]', yj[itr], intercept, cj0)
            g = yj - Yo' * wj
            g = g :- cj0
            Pi[j] = mean(g[ipo]) - mean(g[iva])
            b = qrsolve(XQ, g[iva] \ g[ipo])
            Cj[j] = 2 * abs(b[3])
        }
        Cdon = max(Cj)
        if (dbg) printf("[7] placebos done, Cdon %f\n", Cdon)
        st_numscalar("r(Cdon)", Cdon); st_numscalar("r(above)", sum(Cj :> Cstar))
        st_numscalar("r(lb_Cdon)", theta - Cdon * kap - tq * se); st_numscalar("r(ub_Cdon)", theta + Cdon * kap + tq * se)
        st_numscalar("r(Pimax)", max(abs(Pi)))
        st_matrix("r(placebo)", (Pi, Cj))
        st_matrixcolstripe("r(placebo)", (J(2, 1, ""), ("Pi" \ "curvature")))
    }

    // breakdown persistence on the grid of phibar values
    phs = sort(st_matrix(phm)', 1)
    phistar = .
    knv = J(rows(phs), 1, .)
    if (abs(theta) - tq * se <= 0) phistar = 0
    else {
        for (k = 1; k <= rows(phs); k++) {
            kn = _ssc_knum(lam, xP, xQ, phs[k], alpha, R)
            knv[k] = kn
            if (abs(theta) - kn * se <= 0) {
                phistar = phs[k]
                break
            }
        }
    }
    st_numscalar("r(phistar)", phistar); st_numscalar("r(phimax)", max(phs))
    if (dbg) printf("[8] persistence done\n")

    // Theorem 3: placebo-calibrated drift bound with a simultaneous critical value
    if (doexch) {
        hU = mean(xQ) - mean(xP); G = invsym(A' * A)
        hs = rangen(0, hU, 201); Lm = J(201, n, .)
        for (k = 1; k <= 201; k++) Lm[k, .] = (A * G * (1 \ (mean(xQ) - hs[k])))'
        aa = sqrt(1 / m :+ rowsum(Lm :^ 2))
        if (dbg) printf("[9] family %f x %f, aa %f x %f\n", rows(Lm), cols(Lm), rows(aa), cols(aa))
        Hm = I(n) - A * pinv(A)
        sup = J(0, 1, .)
        for (r0 = 1; r0 <= R; r0 = r0 + 10000) {
            nr = min((10000, R - r0 + 1))
            EP = rnormal(nr, n, 0, 1); EQ = rnormal(nr, 1, 0, 1) / sqrt(m)
            OM = sqrt(rowsum((EP * Hm') :^ 2) / (n - 2))
            sup = sup \ rowmax(abs(EQ :- EP * Lm') :/ (OM * aa'))
        }
        if (dbg) printf("[10] sup %f x %f\n", rows(sup), cols(sup))
        sup = sort(sup, 1); ksim = sup[ceil((1 - alpha) * R)]
        res2 = uP - A * pinv(A) * uP; om2 = sqrt(res2' * res2 / (n - 2))
        th = mean(uQ) :- Lm * uP
        Mhat = (max(abs(Pi)) + ksim * om2 * aa[201]) / abs(hU)
        if (dbg) printf("[11] hs %f x %f, th %f x %f\n", rows(hs), cols(hs), rows(th), cols(th))
        half = Mhat * abs(hs) + ksim * om2 * aa
        ii = .; ww = .
        minindex(half, 1, ii, ww); i = ii[1]
        st_numscalar("r(ksim)", ksim); st_numscalar("r(Mhat)", Mhat); st_numscalar("r(exch_h)", hs[i])
        st_numscalar("r(exch_theta)", th[i]); st_numscalar("r(exch_lb)", th[i] - half[i]); st_numscalar("r(exch_ub)", th[i] + half[i])
        st_numscalar("r(exch_guarantee)", 1 - alpha - 1 / (nd + 1))
    }

    // the whole persistence grid, for the figure. Drawn after Theorem 3 so that
    // every other result is identical with and without the figure
    if (dofull) {
        for (k = 1; k <= rows(phs); k++) if (knv[k] >= .) knv[k] = _ssc_knum(lam, xP, xQ, phs[k], alpha, R)
        PS = ((0 \ phs), (tq \ knv))
        PS = (PS, theta :- PS[., 2] * se, theta :+ PS[., 2] * se)
        st_matrix("r(persist)", PS)
        st_matrixcolstripe("r(persist)", (J(4, 1, ""), ("phibar" \ "crit" \ "lb" \ "ub")))
    }
}
end
