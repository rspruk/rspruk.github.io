*! easc 1.0.0  8oct2026  Rok Spruk
*! Exposure-augmented synthetic control under out-of-span common shocks
*! Implements the diagnostics, estimators and inference of
*! "Synthetic Control under Out-of-Span Common Shocks: Diagnosis, Exposure Balance and Correction"

program define easc, eclass
    version 16.0
    syntax varname(numeric) [if] [in], ID(varname numeric) Time(varname numeric) ///
        TREated(numlist max=1) TRPeriod(numlist max=1) EXPosure(varname numeric)   ///
        [SHock(varname numeric) PHases(numlist) Factors(integer 2) CONtrols(integer 3) ///
         NOCONtrols PLACebo Mlevels(numlist) FRontier(integer 21) Level(cilevel) ///
         SAVing(string) replace GRaph(string) GRAPHPrefix(string) SCHeme(string)]

    marksample touse
    markout `touse' `id' `time' `exposure'
    if "`shock'" != "" markout `touse' `shock'
    local y `varlist'

    * ---------------- checks
    qui count if `touse' & `id' == `treated'
    if r(N) == 0 {
        di as err "treated unit `treated' not found in `id'"
        exit 198
    }
    if "`phases'" == "" {
        qui su `time' if `touse'
        local phases `trperiod' `r(max)'
    }
    local nph : word count `phases'
    if mod(`nph', 2) != 0 {
        di as err "phases() must list start-end pairs, e.g. phases(2012 2015 2016 2017 2018 2024)"
        exit 198
    }
    if "`mlevels'" == "" local mlevels 0 0.5 1
    if "`nocontrols'" != "" local controls 0

    * ---------------- panel checks: balanced, exposure time-invariant
    tempvar sdexp
    qui bysort `id' (`time'): gen double `sdexp' = `exposure' - `exposure'[1] if `touse'
    qui su `sdexp' if `touse'
    if r(max) > 1e-10 | r(min) < -1e-10 {
        di as err "exposure() must be constant within unit (a pre-determined index)"
        exit 198
    }

    * ---------------- run
    tempname W PHI SEPHI GAP GEA DIAG PHT INF FRT LEV PGAP SENS
    mata: easc_main("`y'", "`id'", "`time'", "`exposure'", "`shock'", "`touse'", ///
        `treated', `trperiod', "`phases'", `factors', `controls', ///
        "`placebo'" != "", "`mlevels'", `frontier', `level', ///
        "`W'", "`PHI'", "`SEPHI'", "`GAP'", "`GEA'", "`DIAG'", "`PHT'", "`INF'", "`FRT'", ///
        "`LEV'", "`PGAP'", "`SENS'")

    * ---------------- post
    tempname b
    matrix `b' = `W'
    ereturn post `b', depname(`y') esample(`touse')
    ereturn matrix w        = `W', copy
    ereturn matrix phi      = `PHI', copy
    ereturn matrix se_phi   = `SEPHI', copy
    ereturn matrix gap_sc   = `GAP', copy
    ereturn matrix gap_easc = `GEA', copy
    ereturn matrix diag     = `DIAG', copy
    ereturn matrix phases   = `PHT', copy
    ereturn matrix frontier = `FRT', copy
    if "`placebo'" != "" {
        ereturn matrix inference = `INF', copy
        ereturn matrix placebo_gaps = `PGAP', copy
        ereturn matrix sensitivity = `SENS', copy
    }
    if "`shock'" != "" ereturn matrix leverage = `LEV', copy
    ereturn local cmd "easc"
    ereturn local depvar "`y'"
    ereturn local exposure "`exposure'"
    ereturn scalar treated  = `treated'
    ereturn scalar trperiod = `trperiod'
    ereturn scalar factors  = `factors'

    if "`saving'" != "" {
        preserve
        drop _all
        mata: easc_save("`GAP'", "`GEA'")
        qui save "`saving'", `replace'
        restore
        di as txt "annual gaps saved to `saving'"
    }

    * ---------------- graphs
    if "`graph'" != "" {
        if "`graph'" == "all" local graph gaps leverage frontier sensitivity
        local nph = `nph' / 2
        foreach g of local graph {
            if !inlist("`g'", "gaps", "leverage", "frontier", "sensitivity") {
                di as err "graph(): unknown graph `g'; choose from gaps leverage frontier sensitivity all"
                exit 198
            }
            if "`g'" == "leverage" & "`shock'" == "" {
                di as txt "graph(leverage) needs shock(); skipped"
                continue
            }
            if "`g'" == "sensitivity" & ("`placebo'" == "" | "`shock'" == "") {
                di as txt "graph(sensitivity) needs placebo and shock(); skipped"
                continue
            }
            preserve
            drop _all
            easc_graph_`g', gap(`GAP') gea(`GEA') pgap(`PGAP') lev(`LEV') frt(`FRT') sens(`SENS') pht(`PHT') ///
                trperiod(`trperiod') placebo(`placebo') name(`graphprefix'`g') nph(`nph') scheme(`scheme')
            restore
        }
    }
end


* ---------------- graph helpers (each receives matrices and builds a temporary dataset)
program define easc_graph_gaps
    syntax, name(string) [gap(name) gea(name) pgap(name) lev(name) frt(name) sens(name) pht(name) trperiod(real 0) placebo(string) nph(integer 1) scheme(string)]
    if "`scheme'" != "" local options scheme(`scheme')
    else local options
    qui svmat double `gap', names(g)
    rename (g1 g2) (year gap_sc)
    tempname E
    matrix `E' = `gea'[1..., 2]
    qui svmat double `E', names(e)
    rename e1 gap_easc
    local plac
    local J 0
    if "`placebo'" != "" {
        qui svmat double `pgap', names(p)
        local J = colsof(`pgap')
        forvalues j = 1/`J' {
            local plac `plac' (line p`j' year, lcolor(gs13) lwidth(vthin))
        }
    }
    twoway `plac' (line gap_sc year, lpattern(dash) lcolor(gs5) lwidth(medthick)) ///
        (line gap_easc year, lcolor(navy) lwidth(medthick)), ///
        yline(0, lcolor(gs8)) xline(`=`trperiod'-0.5', lcolor(gs8) lpattern(dot)) ///
        legend(order(`=`J'+1' "Synthetic control" `=`J'+2' "Exposure-augmented SC") rows(1) position(6)) ///
        ytitle("Gap, outcome units") xtitle("") title("Treated unit minus synthetic unit", size(medium)) ///
        name(`name', replace) `options'
end

program define easc_graph_leverage
    syntax, name(string) [gap(name) gea(name) pgap(name) lev(name) frt(name) sens(name) pht(name) trperiod(real 0) placebo(string) nph(integer 1) scheme(string)]
    if "`scheme'" != "" local options scheme(`scheme')
    else local options
    qui svmat double `lev', names(l)
    rename (l1 l2 l3 l4) (year ell kappa kappa_lat)
    qui gen byte post = year >= `trperiod'
    twoway (bar ell year if !post, color(gs9) barwidth(0.8)) (bar ell year if post, color(navy) barwidth(0.8)), ///
        yline(2 -2, lcolor(gs6) lpattern(dash)) yline(0, lcolor(gs8)) ///
        legend(order(1 "pre-treatment" 2 "post-treatment") rows(1) position(6)) ///
        ytitle("Partial leverage of the shock") xtitle("") ///
        title("Silent-factor diagnostic: |values| above 2 after treatment", size(medium)) ///
        name(`name'_ell, replace) `options'
    twoway (line kappa year, lcolor(black) lwidth(medthick)) (line kappa_lat year, lcolor(navy) lpattern(dash)), ///
        xline(`=`trperiod'-0.5', lcolor(gs8) lpattern(dot)) ///
        legend(order(1 "all factors incl. shock" 2 "latent factors only") rows(1) position(6)) ///
        ytitle("Leverage {&kappa}{sub:t}") xtitle("") title("Leverage of the common factors", size(medium)) ///
        name(`name'_kappa, replace) `options'
end

program define easc_graph_frontier
    syntax, name(string) [gap(name) gea(name) pgap(name) lev(name) frt(name) sens(name) pht(name) trperiod(real 0) placebo(string) nph(integer 1) scheme(string)]
    if "`scheme'" != "" local options scheme(`scheme')
    else local options
    qui svmat double `frt', names(f)
    rename (f1 f2 f3) (ei rmspe position)
    qui gen double absei = abs(ei)
    sort absei
    twoway (connected rmspe absei, lcolor(navy) mcolor(navy) msize(small)), ///
        ytitle("Pre-treatment RMSPE") xtitle("Exposure imbalance |EI(w)|") ///
        title("Exposure-fit frontier (right end = SC, left end = EBSC-min)", size(medium)) ///
        name(`name', replace) `options'
end

program define easc_graph_sensitivity
    syntax, name(string) [gap(name) gea(name) pgap(name) lev(name) frt(name) sens(name) pht(name) trperiod(real 0) placebo(string) nph(integer 1) scheme(string)]
    if "`scheme'" != "" local options scheme(`scheme')
    else local options
    qui svmat double `sens', names(s)
    rename (s1 s2 s3 s4 s5 s6) (est phase M lo hi tau)
    forvalues p = 1/`nph' {
        twoway (rarea lo hi M if est == 1 & phase == `p', color(gs10%50) lwidth(none)) ///
               (rarea lo hi M if est == 2 & phase == `p', color(navy%40) lwidth(none)) ///
               (line tau M if est == 1 & phase == `p', lcolor(gs5) lpattern(dash)) ///
               (line tau M if est == 2 & phase == `p', lcolor(navy)), ///
            yline(0, lcolor(gs8)) xtitle("M (share of pre-treatment misfit treated as systematic)") ///
            ytitle("Bias-aware confidence interval") ///
            legend(order(1 "SC" 2 "EASC") rows(1) position(6)) ///
            title("Phase `p': sensitivity to systematic misfit", size(medium)) ///
            name(`name'`p', replace) `options'
    }
end


* =====================================================================================
version 16.0
mata:
mata clear
mata set matastrict on

// ---------------------------------------------------------------- simplex least squares
// min ||a - A w||^2 s.t. w >= 0, sum w = 1 [, s0'w = starget]
// Lawson-Hanson NNLS on an augmented system (heavy rows impose the equalities)
real colvector easc_nnls(real colvector b, real matrix A, real scalar maxit)
{
    real scalar n, m, i, j, alpha, t, it
    real colvector x, w, z, P
    real rowvector pidx
    real matrix Ap

    m = rows(A)
    n = cols(A)
    x = J(n, 1, 0)
    P = J(n, 1, 0)
    it = 0
    while (it < maxit) {
        it = it + 1
        w = A' * (b - A * x)
        if (sum(P) == n) {
            break
        }
        j = 0
        t = 1e-10
        for (i = 1; i <= n; i++) {
            if (P[i] == 0) {
                if (w[i] > t) {
                    t = w[i]
                    j = i
                }
            }
        }
        if (j == 0) {
            break
        }
        P[j] = 1
        while (1) {
            pidx = select((1..n), P')
            Ap = A[., pidx]
            z = J(n, 1, 0)
            z[pidx] = qrsolve(Ap, b)
            if (min(z[pidx]) > 0) {
                x = z
                break
            }
            alpha = 1
            for (i = 1; i <= n; i++) {
                if (P[i] == 1) {
                    if (z[i] <= 0) {
                        if (x[i] / (x[i] - z[i]) < alpha) {
                            alpha = x[i] / (x[i] - z[i])
                        }
                    }
                }
            }
            x = x + alpha * (z - x)
            for (i = 1; i <= n; i++) {
                if (P[i] == 1) {
                    if (x[i] <= 1e-12) {
                        P[i] = 0
                        x[i] = 0
                    }
                }
            }
        }
    }
    return(x)
}

real colvector easc_simplex(real colvector a, real matrix A, | real colvector s0, real scalar starget)
{
    real scalar M, nJ
    real matrix Aa
    real colvector ba, w
    nJ = cols(A)
    M = 1e3 * (max(abs(A)) + 1)
    Aa = A \ (M * J(1, nJ, 1))
    ba = a \ M
    if (args() > 2) {
        Aa = Aa \ (M * s0')
        ba = ba \ (M * starget)
    }
    w = easc_nnls(ba, Aa, 50 * nJ)
    w = w :* (w :> 0)
    return(w / sum(w))
}

// ---------------------------------------------------------------- latent factors from donors
real matrix easc_factors(real matrix Y0, real scalar k, real colvector s0)
{
    real matrix Z, U, Vt
    real colvector sv, sc
    Z = Y0 :- mean(Y0)                       // remove unit means (rows = periods)
    Z = Z :- mean(Z')'                       // remove period means
    sc = (s0 :- mean(s0)) / sqrt(variance(s0))
    Z = Z - (Z * sc / (sc' * sc)) * sc'      // remove cross-sectional projection on exposure, period by period
    svd(Z, U = ., sv = ., Vt = .)
    return(U[., 1..k] :* (sv[1..k]' / sqrt(rows(Z))))
}

// leverage of each row of G relative to the pre-period second-moment matrix
real colvector easc_leverage(real matrix G, real colvector pre)
{
    real matrix S0, Si
    S0 = (select(G, pre))' * select(G, pre) / sum(pre)
    Si = pinv(S0)
    return(sqrt(rowsum((G * Si) :* G)))
}

real scalar easc_phase_leverage(real matrix G, real colvector pre, real colvector mask)
{
    real matrix S0
    real rowvector gb
    S0 = (select(G, pre))' * select(G, pre) / sum(pre)
    gb = mean(select(G, mask))
    return(sqrt(gb * pinv(S0) * gb'))
}

// partial leverage of the observed shock x on factors F (fitted on the pre-period)
real colvector easc_partial(real colvector x, real matrix F, real colvector pre, real scalar r2)
{
    real matrix X
    real colvector b, res
    X = J(rows(x), 1, 1), F
    b = qrsolve(select(X, pre), select(x, pre))
    res = x - X * b
    r2 = 1 - variance(select(res, pre)) / variance(select(x, pre))
    return(res / sqrt(mean(select(res, pre) :^ 2)))
}

// ---------------------------------------------------------------- exposure regression
// phi_t: coefficient on s_j in the cross-donor regression of Y_jt - Y_jb on (1, s_j, c_j);
// c_j = growth within each of `nctl` sub-periods of the pre-treatment years before b
void easc_phi(real matrix Y0, real colvector s0, real colvector pre, real scalar bidx, real scalar nctl,
              real colvector phi, real colvector se)
{
    real scalar T, nJ, j, k, n, lo, hi
    real matrix D, C, X, JK, b
    real colvector preidx
    T = rows(Y0); nJ = cols(Y0)
    D = Y0 :- Y0[bidx, .]                                // T x nJ
    preidx = select((1::T), pre :& ((1::T) :< bidx))
    n = rows(preidx)
    C = J(nJ, 0, .)
    hi = 0
    for (k = 1; k <= nctl; k++) {                         // same split as numpy.array_split
        lo = hi + 1
        hi = hi + floor(n / nctl) + (k <= mod(n, nctl))
        C = C, (Y0[preidx[hi], .] - Y0[preidx[lo], .])'
    }
    X = J(nJ, 1, 1), s0, C
    b = qrsolve(X, D')                                    // (2+nctl) x T
    phi = b[2, .]'
    JK = J(nJ, T, .)
    for (j = 1; j <= nJ; j++) {
        b = qrsolve(X[select(1::nJ, (1::nJ) :!= j), .], (D[., select(1::nJ, (1::nJ) :!= j)])')
        JK[j, .] = b[2, .]
    }
    se = sqrt((nJ - 1) / nJ * colsum((JK :- mean(JK)) :^ 2))'
}

// ---------------------------------------------------------------- bias-aware critical value
real scalar easc_cv(real scalar b, real scalar alpha)
{
    if (abs(b) < 1e-10) {
        return(invnormal(1 - alpha / 2))
    }
    return(sqrt(invnchi2(1, b^2, 1 - alpha)))
}

real scalar easc_breakdown(real scalar tau, real scalar kbar, real scalar rmspe, real scalar sd, real scalar alpha)
{
    real scalar lo, hi, mid, i
    if (abs(tau) <= easc_cv(0, alpha) * sd) {
        return(0)
    }
    lo = 0
    hi = 5
    if (abs(tau) > easc_cv(hi * kbar * rmspe / sd, alpha) * sd) {
        return(hi)
    }
    for (i = 1; i <= 60; i++) {
        mid = (lo + hi) / 2
        if (abs(tau) > easc_cv(mid * kbar * rmspe / sd, alpha) * sd) {
            lo = mid
        }
        else {
            hi = mid
        }
    }
    return(lo)
}

// ---------------------------------------------------------------- one unit (treated or placebo)
struct easc_unit {
    real colvector w, g, ge, phi, se
    real scalar ei, rmspe
    real matrix Fl
}

struct easc_unit scalar easc_fit(real colvector y1, real matrix Y0, real scalar s1, real colvector s0,
                                 real colvector pre, real scalar bidx, real scalar k, real scalar nctl,
                                 real colvector x)
{
    struct easc_unit scalar u
    real colvector phi, se
    u.w = easc_simplex(select(y1, pre), select(Y0, pre))
    u.g = y1 - Y0 * u.w
    u.ei = s1 - s0' * u.w
    if (nctl >= 0) {
        easc_phi(Y0, s0, pre, bidx, nctl, phi = ., se = .)
        u.phi = phi; u.se = se
        u.ge = u.g - phi * u.ei
    }
    u.Fl = easc_factors(Y0, k, s0)
    u.rmspe = sqrt(mean(select(u.g, pre) :^ 2))
    return(u)
}

// ---------------------------------------------------------------- main
void easc_main(string scalar y, string scalar id, string scalar time, string scalar sexp, string scalar sshock,
               string scalar touse, real scalar treated, real scalar trp, string scalar sphases,
               real scalar k, real scalar nctl, real scalar doplac, string scalar smlev, real scalar nfr,
               real scalar level,
               string scalar Wn, string scalar PHIn, string scalar SEn, string scalar GAPn, string scalar GEAn,
               string scalar DIAGn, string scalar PHTn, string scalar INFn, string scalar FRTn,
               string scalar LEVn, string scalar PGAPn, string scalar SENSn)
{
    real matrix D, Y, Y0, G, Gl, FRT, PHT, INF, DIAG, PGAP, SENS
    real scalar Mg, nMg
    real colvector yrs, uid, s, y1, s0, pre, x, phases, mask, ell, z, mlev
    real scalar nJ, T, s1, bidx, i, j, p, r2, dmin, alpha, tau, sd, kb, Mstar, ext, hw, m, pval, a, b
    real rowvector tmp
    real colvector val
    string scalar lab
    struct easc_unit scalar u
    struct easc_unit colvector plac

    // ----- data to wide T x (nJ+1), treated first
    if (sshock != "") {
        D = st_data(., (id, time, y, sexp, sshock), touse)
    }
    else {
        D = st_data(., (id, time, y, sexp), touse)
    }
    uid = uniqrows(D[., 1]); yrs = uniqrows(D[., 2])
    uid = treated \ select(uid, uid :!= treated)
    nJ = rows(uid) - 1; T = rows(yrs)
    Y = J(T, nJ + 1, .); s = J(nJ + 1, 1, .)
    for (j = 1; j <= nJ + 1; j++) {
        for (i = 1; i <= T; i++) {
            val = select(D[., 3], (D[., 1] :== uid[j]) :& (D[., 2] :== yrs[i]))
            if (rows(val) != 1) {
                _error(3000, "panel must be balanced: one observation per unit and period")
            }
            Y[i, j] = val[1]
        }
        val = select(D[., 4], D[., 1] :== uid[j])
        s[j] = val[1]
    }
    x = J(T, 1, .)
    if (sshock != "") {
        for (i = 1; i <= T; i++) {
            val = select(D[., 5], (D[., 1] :== treated) :& (D[., 2] :== yrs[i]))
            x[i] = val[1]
        }
    }
    y1 = Y[., 1]; Y0 = Y[., 2..nJ+1]; s1 = s[1]; s0 = s[2..nJ+1]
    pre = yrs :< trp
    bidx = max(select(1::T, pre))
    phases = strtoreal(tokens(sphases))'
    mlev = strtoreal(tokens(smlev))'
    alpha = 1 - level / 100

    // ----- treated unit
    u = easc_fit(y1, Y0, s1, s0, pre, bidx, k, nctl, x)
    Gl = J(T, 1, 1), u.Fl
    dmin = max((0, s1 - max(s0), min(s0) - s1))

    // ----- frontier: exposure targets from SC imbalance down to the attainable minimum
    FRT = J(nfr, 3, .)
    for (i = 1; i <= nfr; i++) {
        p = (i - 1) / (nfr - 1)                                    // 0 = SC end, 1 = balanced end
        a = s0' * u.w; b = s1 - sign(u.ei) * dmin                  // exposure of SC, and of EBSC-min
        tmp = easc_simplex(select(y1, pre), select(Y0, pre), s0, a + p * (b - a))'
        FRT[i, .] = (s1 - s0' * tmp', sqrt(mean(select(y1 - Y0 * tmp', pre) :^ 2)), p)
    }

    // ----- diagnostics
    DIAG = J(10, 1, .)
    DIAG[1] = s1; DIAG[2] = min(s0); DIAG[3] = max(s0); DIAG[4] = dmin; DIAG[5] = u.ei
    DIAG[6] = u.rmspe; DIAG[7] = FRT[nfr, 2]
    if (sshock != "") {
        ell = easc_partial(x, u.Fl, pre, r2 = .)
        G = Gl, x
        DIAG[8] = r2; DIAG[9] = max(abs(select(ell, pre))); DIAG[10] = max(abs(select(ell, !pre)))
        st_matrix(LEVn, (yrs, ell, easc_leverage(G, pre), easc_leverage(Gl, pre)))
    }

    // ----- phase table: SC, EASC, EBSC-min, plus leverage
    PHT = J(rows(phases) / 2, 8, .)
    for (p = 1; p <= rows(phases) / 2; p++) {
        mask = (yrs :>= phases[2*p-1]) :& (yrs :<= phases[2*p])
        PHT[p, 1] = phases[2*p-1]; PHT[p, 2] = phases[2*p]
        PHT[p, 3] = mean(select(u.g, mask))
        PHT[p, 4] = mean(select(u.ge, mask))
        tmp = easc_simplex(select(y1, pre), select(Y0, pre), s0, s1 - sign(u.ei) * dmin)'
        PHT[p, 5] = mean(select(y1 - Y0 * tmp', mask))
        PHT[p, 6] = mean(select(u.phi, mask))
        PHT[p, 7] = (sshock != "" ? easc_phase_leverage(G, pre, mask) : .)
        PHT[p, 8] = (sshock != "" ? sqrt(easc_phase_leverage(Gl, pre, mask)^2 + ell[bidx]^2) : .)
    }

    // ----- placebo inference (SC and EASC): standardised placebo dispersion, bias-aware CIs, breakdown
    INF = J(0, 8 + 2 * rows(mlev), .)
    SENS = J(0, 6, .)
    nMg = 51
    if (doplac) {
        plac = easc_unit(nJ, 1)
        PGAP = J(T, nJ, .)
        for (j = 1; j <= nJ; j++) {
            plac[j] = easc_fit(Y0[., j], Y0[., select(1::nJ, (1::nJ) :!= j)], s0[j],
                               s0[select(1::nJ, (1::nJ) :!= j)], pre, bidx, k, nctl, x)
            PGAP[., j] = plac[j].ge
        }
        st_matrix(PGAPn, PGAP)
        for (m = 1; m <= 2; m++) {                                    // 1 = SC, 2 = EASC
            for (p = 1; p <= rows(phases) / 2; p++) {
                mask = (yrs :>= phases[2*p-1]) :& (yrs :<= phases[2*p])
                tau = (m == 1 ? PHT[p, 3] : PHT[p, 4])
                z = J(nJ, 1, .); pval = 1
                for (j = 1; j <= nJ; j++) {
                    z[j] = mean(select((m == 1 ? plac[j].g : plac[j].ge), mask)) / plac[j].rmspe
                    if (abs(z[j]) >= abs(tau) / u.rmspe) {
                        pval = pval + 1
                    }
                }
                pval = pval / (nJ + 1)
                sd = sqrt(variance(z)) * u.rmspe
                if (m == 2) {
                    ext = abs(u.ei) * mean(select(u.se, mask))
                    sd = sqrt(sd^2 + ext^2)
                }
                kb = (sshock != "" ? (m == 1 ? PHT[p, 7] : PHT[p, 8]) : .)
                tmp = (m, phases[2*p-1], phases[2*p], tau, pval, sd, kb)
                for (i = 1; i <= rows(mlev); i++) {
                    hw = (kb < . ? easc_cv(mlev[i] * kb * u.rmspe / sd, alpha) * sd : easc_cv(0, alpha) * sd)
                    tmp = tmp, tau - hw, tau + hw
                }
                Mstar = (kb < . ? easc_breakdown(tau, kb, u.rmspe, sd, alpha) : .)
                INF = INF \ (tmp, Mstar)
                if (kb < .) {
                    for (i = 1; i <= nMg; i++) {
                        Mg = (i - 1) / (nMg - 1)
                        hw = easc_cv(Mg * kb * u.rmspe / sd, alpha) * sd
                        SENS = SENS \ (m, p, Mg, tau - hw, tau + hw, tau)
                    }
                }
            }
        }
    }

    // ----- display
    printf("\n{txt}Exposure-augmented synthetic control{res}\n")
    printf("{txt}Treated unit %g, first treated period %g, %g donors, %g pre-treatment periods\n", treated, trp, nJ, sum(pre))
    printf("\n{txt}Exposure geometry\n")
    printf("{txt}  treated exposure s1          = {res}%9.3f\n", s1)
    printf("{txt}  donor exposure range         = {res}[%8.3f, %8.3f]\n", min(s0), max(s0))
    printf("{txt}  distance to hull d_min       = {res}%9.3f\n", dmin)
    printf("{txt}  exposure imbalance EI(w_SC)  = {res}%9.3f\n", u.ei)
    printf("{txt}  pre-treatment RMSPE: SC      = {res}%9.4f{txt}   EBSC-min = {res}%9.4f\n", u.rmspe, FRT[nfr, 2])
    if (sshock != "") {
        printf("\n{txt}Silent-factor diagnostic for %s (%g latent factors)\n", sshock, k)
        printf("{txt}  pre-treatment R2 of shock on factors   = {res}%6.3f\n", r2)
        printf("{txt}  max |partial leverage|, pre / post     = {res}%6.2f{txt} / {res}%6.2f\n", DIAG[9], DIAG[10])
    }
    printf("\n{txt}{hline 78}\n")
    printf("{txt}%-14s %10s %10s %10s %10s %8s %8s\n", "Phase", "SC", "EASC", "EBSC-min", "100*phi", "kappa", "kappaEA")
    printf("{txt}{hline 78}\n")
    for (p = 1; p <= rows(PHT); p++) {
        printf("{res}%5.0f-%-8.0f %10.4f %10.4f %10.4f %10.3f %8.2f %8.2f\n", PHT[p,1], PHT[p,2], PHT[p,3], PHT[p,4], PHT[p,5], 100*PHT[p,6], PHT[p,7], PHT[p,8])
    }
    printf("{txt}{hline 78}\n")
    if (doplac) {
        printf("\n{txt}Placebo inference and bias-aware intervals (%g%% level; bias bound M*kappa*RMSPE)\n", level)
        printf("{txt}{hline 78}\n")
        printf("{txt}%-5s %-11s %9s %7s %8s", "Est.", "Phase", "estimate", "p", "noise sd")
        for (i = 1; i <= rows(mlev); i++) {
            printf(" %17s", "CI M=" + strofreal(mlev[i]))
        }
        printf(" %6s\n", "M*")
        printf("{txt}{hline 78}\n")
        for (i = 1; i <= rows(INF); i++) {
            lab = (INF[i,1] == 1 ? "SC" : "EASC")
            printf("{res}%-5s %5.0f-%-5.0f %9.4f %7.3f %8.4f", lab, INF[i,2], INF[i,3], INF[i,4], INF[i,5], INF[i,6])
            for (m = 1; m <= rows(mlev); m++) {
                printf(" [%7.3f,%7.3f]", INF[i, 6+2*m], INF[i, 7+2*m])
            }
            printf(" %6.2f\n", INF[i, cols(INF)])
        }
        printf("{txt}{hline 78}\n")
    }

    // ----- return
    st_matrix(Wn, u.w'); st_matrixcolstripe(Wn, (J(nJ, 1, ""), "d" :+ strofreal(uid[2..nJ+1])))
    st_matrix(PHIn, (yrs, u.phi)); st_matrix(SEn, (yrs, u.se))
    st_matrix(GAPn, (yrs, u.g)); st_matrix(GEAn, (yrs, u.ge))
    st_matrix(DIAGn, DIAG)
    st_matrixrowstripe(DIAGn, (J(10, 1, ""), ("s1" \ "s0_min" \ "s0_max" \ "d_min" \ "EI_sc" \ "rmspe_sc" \ "rmspe_ebscmin" \ "r2_shock" \ "maxabs_ell_pre" \ "maxabs_ell_post")))
    st_matrix(PHTn, PHT)
    st_matrixcolstripe(PHTn, (J(8, 1, ""), ("start" \ "end" \ "sc" \ "easc" \ "ebsc_min" \ "phi" \ "kappa_sc" \ "kappa_easc")))
    st_matrix(FRTn, FRT); st_matrixcolstripe(FRTn, (J(3, 1, ""), ("EI" \ "rmspe" \ "position")))
    if (doplac) {
        st_matrix(INFn, INF)
        st_matrix(SENSn, SENS)
    }
}

void easc_save(string scalar GAPn, string scalar GEAn)
{
    real matrix G, E
    G = st_matrix(GAPn); E = st_matrix(GEAn)
    st_addobs(rows(G))
    (void) st_addvar("double", ("year", "gap_sc", "gap_easc"))
    st_store(., ("year", "gap_sc", "gap_easc"), (G, E[., 2]))
}

end
