* splitsc_examples.do  The three classic applications of the paper.
* Run from the folder that contains splitsc.ado and the data subfolder.
clear all
set more off
adopath ++ "`c(pwd)'"

* California, Proposition 99 (Abadie, Diamond and Hainmueller 2010)
import delimited using "prop99.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("California") trperiod(1989) trainend(1979) seed(2026)

* West Germany, reunification (Abadie, Diamond and Hainmueller 2015), with Theorem 3
import delimited using "germany.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exch seed(2026)
* anticipation check, 1990 treated as an excluded period
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exclude(1990) seed(2026)

* Basque Country (Abadie and Gardeazabal 2003), a stated curvature bound and simplex weights
import delimited using "basque.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("Basque Country (Pais Vasco)") trperiod(1970) trainend(1961) curvature(0.003) nointercept seed(2026)
matrix list r(weights)
