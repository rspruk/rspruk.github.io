* splitsc_cert.do  Certification: compares splitsc with the Python reference implementation used in the paper.
* Run from the folder that contains splitsc.ado and the data subfolder:  do splitsc_cert.do
clear all
set more off
adopath ++ "`c(pwd)'"
local fail 0
program define chk
    args name got want tol
    local d = abs(`got' - `want') / max(1, abs(`want'))
    if `d' > `tol' {
        di as err "FAIL `name': got " %14.8g `got' " expected " %14.8g `want'
        c_local fail = 1
    }
    else di as txt "ok   `name'"
end

* ---- prop99_sc
import delimited using "prop99.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("California") trperiod(1989) trainend(1979) nointercept phibar(0.5) reps(2000) seed(1)
chk prop99_sc_theta r(theta) -5.74687770063732 1e-05
chk prop99_sc_lb r(lb) -12.794196239189375 1e-05
chk prop99_sc_ub r(ub) 1.300440837914734 1e-05
chk prop99_sc_nu r(nu) 7.0 0
chk prop99_sc_kappa r(kappa) 57.749999998735575 1e-06
chk prop99_sc_Cstar r(Cstar) 0.0 1e-05
chk prop99_sc_Cdon r(Cdon) 0.7034722018206057 1e-05
chk prop99_sc_lb_Cdon r(lb_Cdon) -53.419715893439864 1e-05
chk prop99_sc_ub_Cdon r(ub_Cdon) 41.925960492165224 1e-05
chk prop99_sc_Pimax r(Pimax) 20.422608165882924 1e-05
chk prop99_sc_above r(above) 38.0 0

* ---- prop99_dscm
import delimited using "prop99.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("California") trperiod(1989) trainend(1979) phibar(0.5) reps(2000) seed(1)
chk prop99_dscm_theta r(theta) -8.054470317119168 1e-05
chk prop99_dscm_lb r(lb) -15.137030068054194 1e-05
chk prop99_dscm_ub r(ub) -0.9719105661841416 1e-05
chk prop99_dscm_nu r(nu) 7.0 0
chk prop99_dscm_kappa r(kappa) 57.749999998735575 1e-06
chk prop99_dscm_Cstar r(Cstar) 0.016829620194033272 1e-05
chk prop99_dscm_Cdon r(Cdon) 0.7001531776909897 1e-05
chk prop99_dscm_lb_Cdon r(lb_Cdon) -55.57087607882356 1e-05
chk prop99_dscm_ub_Cdon r(ub_Cdon) 39.46193544458522 1e-05
chk prop99_dscm_Pimax r(Pimax) 34.38379235909416 1e-05
chk prop99_dscm_above r(above) 35.0 0

* ---- germany_sc
import delimited using "germany.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) nointercept phibar(0.5) reps(2000) seed(1)
chk germany_sc_theta r(theta) -1759.6278355752277 1e-05
chk germany_sc_lb r(lb) -2161.2961686619733 1e-05
chk germany_sc_ub r(ub) -1357.959502488482 1e-05
chk germany_sc_nu r(nu) 13.0 0
chk germany_sc_kappa r(kappa) 95.66666667007985 1e-06
chk germany_sc_Cstar r(Cstar) 14.194698631777348 1e-05
chk germany_sc_Cdon r(Cdon) 30.39145720180204 1e-05
chk germany_sc_lb_Cdon r(lb_Cdon) -5068.745574404767 1e-05
chk germany_sc_ub_Cdon r(ub_Cdon) 1549.4899032543112 1e-05
chk germany_sc_Pimax r(Pimax) 5547.88857141072 1e-05
chk germany_sc_above r(above) 6.0 0

* ---- germany_dscm
import delimited using "germany.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) phibar(0.5) reps(2000) seed(1)
chk germany_dscm_theta r(theta) -2112.113555477692 1e-05
chk germany_dscm_lb r(lb) -2398.1749660765613 1e-05
chk germany_dscm_ub r(ub) -1826.0521448788227 1e-05
chk germany_dscm_nu r(nu) 13.0 0
chk germany_dscm_kappa r(kappa) 95.66666667007985 1e-06
chk germany_dscm_Cstar r(Cstar) 19.087653081676024 1e-05
chk germany_dscm_Cdon r(Cdon) 43.494989684644864 1e-05
chk germany_dscm_lb_Cdon r(lb_Cdon) -6559.195646056043 1e-05
chk germany_dscm_ub_Cdon r(ub_Cdon) 2334.9685351006588 1e-05
chk germany_dscm_Pimax r(Pimax) 5376.504659828561 1e-05
chk germany_dscm_above r(above) 4.0 0

* ---- basque_sc
import delimited using "basque.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("Basque Country (Pais Vasco)") trperiod(1970) trainend(1961) nointercept phibar(0.5) reps(2000) seed(1)
chk basque_sc_theta r(theta) -1.4518553852112563 1e-05
chk basque_sc_lb r(lb) -1.603008452451847 1e-05
chk basque_sc_ub r(ub) -1.3007023179706656 1e-05
chk basque_sc_nu r(nu) 6.0 0
chk basque_sc_kappa r(kappa) 191.99999999836072 1e-06
chk basque_sc_Cstar r(Cstar) 0.00677449123948839 1e-05
chk basque_sc_Cdon r(Cdon) 0.005489093732212166 1e-05
chk basque_sc_lb_Cdon r(lb_Cdon) -2.656914449027585 1e-05
chk basque_sc_ub_Cdon r(ub_Cdon) -0.24679632139492785 1e-05
chk basque_sc_Pimax r(Pimax) 1.2307823557732793 1e-05
chk basque_sc_above r(above) 0.0 0

* ---- basque_dscm
import delimited using "basque.csv", clear varnames(1) asdouble
splitsc y, unit(unit) time(year) treated("Basque Country (Pais Vasco)") trperiod(1970) trainend(1961) phibar(0.5) reps(2000) seed(1)
chk basque_dscm_theta r(theta) -0.8585380378016612 1e-05
chk basque_dscm_lb r(lb) -1.0828055263916168 1e-05
chk basque_dscm_ub r(ub) -0.6342705492117057 1e-05
chk basque_dscm_nu r(nu) 6.0 0
chk basque_dscm_kappa r(kappa) 191.99999999836072 1e-06
chk basque_dscm_Cstar r(Cstar) 0.003303492443839172 1e-05
chk basque_dscm_Cdon r(Cdon) 0.005881405729601288 1e-05
chk basque_dscm_lb_Cdon r(lb_Cdon) -2.2120354264654227 1e-05
chk basque_dscm_ub_Cdon r(ub_Cdon) 0.4949593508621005 1e-05
chk basque_dscm_Pimax r(Pimax) 1.25743405843358 1e-05
chk basque_dscm_above r(above) 5.0 0

if `fail' di as err _n "Some checks failed. Placebo counts can differ by one when a placebo fit has non-unique weights."
else di as txt _n "All checks passed."
