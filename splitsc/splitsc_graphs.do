* splitsc_graphs.do  Figures for the three classic applications.
* Run from the folder that contains splitsc.ado and the data files.
clear all
discard
set more off
set tracedepth 2
adopath ++ "`c(pwd)'"
capture mkdir figures
capture log close _sscg
log using "figures/splitsc_graphs.log", text replace name(_sscg)

* West Germany, all five figures
import delimited using "germany.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exch seed(2026) graph(all) gname(de_)
di "check r() after graphs: Cstar = " r(Cstar) ", Cdon = " r(Cdon) ", lb = " r(lb)
foreach g in paths gap sensitivity placebo persistence {
    graph export "figures/germany_`g'.png", name(de_`g') replace width(1800)
}

* West Germany with 1990 excluded, gap figure only
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exclude(1990) seed(2026) graph(gap) gname(de90_)
graph export "figures/germany1990_gap.png", name(de90_gap) replace width(1800)

* California
import delimited using "prop99.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("California") trperiod(1989) trainend(1979) seed(2026) graph(all) gname(ca_)
foreach g in paths gap sensitivity placebo persistence {
    graph export "figures/california_`g'.png", name(ca_`g') replace width(1800)
}

* Basque Country, a stated curvature bound and no placebos
import delimited using "basque.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("Basque Country (Pais Vasco)") trperiod(1970) trainend(1961) curvature(0.003) nointercept noplacebo seed(2026) graph(paths gap sensitivity)  gname(bq_)
foreach g in paths gap sensitivity {
    graph export "figures/basque_`g'.png", name(bq_`g') replace width(1800)
}
di "all figures done"
log close _sscg
