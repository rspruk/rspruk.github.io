{smcl}
{* *! version 1.0.0  8oct2026}{...}
{title:Title}

{phang}
{bf:easc} {hline 2} Exposure-augmented synthetic control under out-of-span common shocks


{title:Syntax}

{p 8 17 2}
{cmd:easc} {it:depvar} {ifin}{cmd:,} {opt id(varname)} {opt t:ime(varname)} {opt tre:ated(#)} {opt trp:eriod(#)} {opt exp:osure(varname)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt id(varname)}}unit identifier{p_end}
{synopt:{opt t:ime(varname)}}time variable (balanced panel){p_end}
{synopt:{opt tre:ated(#)}}value of {it:id} for the treated unit{p_end}
{synopt:{opt trp:eriod(#)}}first treated period{p_end}
{synopt:{opt exp:osure(varname)}}pre-determined exposure index, constant within unit (e.g. oil rents, % of GDP){p_end}

{syntab:Diagnostics}
{synopt:{opt sh:ock(varname)}}observed common shock (e.g. log oil price), constant across units within period; enables leverage diagnostics and bias-aware intervals{p_end}
{synopt:{opt f:actors(#)}}number of latent factors estimated from the donors; default {cmd:factors(2)}{p_end}

{syntab:Estimation}
{synopt:{opt ph:ases(numlist)}}start-end pairs of post-treatment phases over which effects are averaged; default: whole post-treatment period{p_end}
{synopt:{opt con:trols(#)}}number of pre-treatment growth sub-periods used as controls in the exposure regression; default {cmd:controls(3)}{p_end}
{synopt:{opt nocon:trols}}omit the growth controls{p_end}
{synopt:{opt fr:ontier(#)}}number of points on the exposure-fit frontier; default {cmd:frontier(21)}{p_end}

{syntab:Inference}
{synopt:{opt plac:ebo}}in-space placebos: placebo p-values, noise scale, bias-aware intervals and breakdown values{p_end}
{synopt:{opt m:levels(numlist)}}sensitivity values M for the bias bound M*kappa*RMSPE; default {cmd:mlevels(0 0.5 1)}{p_end}
{synopt:{opt l:evel(#)}}confidence level; default {cmd:level(95)}{p_end}

{syntab:Output}
{synopt:{opt sav:ing(filename)}}save annual SC and EASC gaps to a Stata dataset{p_end}
{synopt:{opt replace}}overwrite {it:filename}{p_end}
{synopt:{opt gr:aph(names)}}draw graphs; any of {cmd:gaps} (SC and EASC gap paths, with placebo gaps if {opt placebo}), {cmd:leverage} (partial leverage of the shock and leverage of the factors; needs {opt shock()}), {cmd:frontier} (exposure-fit frontier), {cmd:sensitivity} (bias-aware intervals as a function of M, one graph per phase; needs {opt placebo} and {opt shock()}), or {cmd:all}{p_end}
{synopt:{opt graphp:refix(string)}}prefix for graph names, e.g. {cmd:graphprefix(iran_)}{p_end}
{synopt:{opt sch:eme(schemename)}}graph scheme{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:easc} implements the diagnostics, estimators and inference of
"Synthetic Control under Out-of-Span Common Shocks: Diagnosis, Exposure Balance and Correction"
(Spruk, 2026). The setting is a synthetic-control design in which an observed common shock
(a commodity price, say) moves with the latent factors before treatment and departs from them
afterwards, so that pre-treatment fit does not pin down the synthetic unit's exposure to the shock.

{pstd}
The command reports, for the treated unit:

{phang2}1. {bf:Exposure geometry}: the treated unit's exposure, the donors' range, the distance to the
donors' exposure hull {it:d_min}, the exposure imbalance EI of the synthetic-control weights, and the
exposure-fit frontier (returned in {cmd:e(frontier)}).

{phang2}2. {bf:Silent-factor diagnostic} (with {opt shock()}): the pre-treatment R-squared of the shock on
the donors' latent factors and the partial leverage of the shock before and after treatment. A high
R-squared with post-treatment |partial leverage| well above 2 indicates a silent factor.

{phang2}3. {bf:Estimates by phase}: synthetic control (SC), the exposure-augmented synthetic control (EASC),
which adds the correction phi_t * EI with phi_t from a cross-donor regression of outcome changes on exposure
and pre-treatment growth, and the exposure-balanced synthetic control at the smallest attainable imbalance
(EBSC-min); the phase average of phi_t; and the leverage kappa of each phase for SC and EASC.

{phang2}4. {bf:Inference} (with {opt placebo}): placebo p-values, the standardised placebo noise scale,
bias-aware confidence intervals for each M in {opt mlevels()}, and the breakdown value M* (the largest M
for which the interval excludes zero).

{pstd}
Weights are simplex least squares solved by non-negative least squares on an augmented system; the
exposure-balanced weights add an exposure equality as a further heavily weighted row. Latent factors are
the leading left singular vectors of two-way demeaned donor outcomes after removing, period by period,
the cross-sectional projection on exposure. Standard errors of phi_t are delete-one-donor jackknife.
With {opt placebo}, every donor is treated in turn with the remaining donors as controls.


{title:Stored results}

{pstd}{cmd:easc} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(w)}}synthetic-control weights (columns named by donor id){p_end}
{synopt:{cmd:e(gap_sc)}}period, SC gap{p_end}
{synopt:{cmd:e(gap_easc)}}period, EASC gap{p_end}
{synopt:{cmd:e(phi)}, {cmd:e(se_phi)}}period, phi_t and its jackknife s.e.{p_end}
{synopt:{cmd:e(phases)}}by phase: start, end, SC, EASC, EBSC-min, mean phi, kappa (SC), kappa (EASC){p_end}
{synopt:{cmd:e(diag)}}exposure geometry and silent-factor diagnostics{p_end}
{synopt:{cmd:e(frontier)}}exposure-fit frontier: EI, RMSPE, position (0 = SC, 1 = EBSC-min){p_end}
{synopt:{cmd:e(inference)}}with {opt placebo}: estimator (1 SC, 2 EASC), phase, estimate, p, noise sd, kappa, CI bounds per M, M*{p_end}
{synopt:{cmd:e(placebo_gaps)}}with {opt placebo}: EASC gaps of every placebo unit (periods by donors){p_end}
{synopt:{cmd:e(sensitivity)}}with {opt placebo} and {opt shock()}: estimator, phase, M, CI lower, CI upper, estimate on a grid of M{p_end}
{synopt:{cmd:e(leverage)}}with {opt shock()}: period, partial leverage of the shock, leverage with and without the shock{p_end}


{title:Example}

{phang}{cmd:. import delimited iran_easc.csv, clear}{p_end}
{phang}{cmd:. gen lbrent = ...}{space 20}(log Brent price by year){p_end}
{phang}{cmd:. easc lgdppc, id(id) time(year) treated(1) trperiod(2012) exposure(oilrent) shock(lbrent) phases(2012 2015 2016 2017 2018 2024) placebo graph(all)}{p_end}


{title:Reference}

{phang}Spruk, R. (2026). Synthetic Control under Out-of-Span Common Shocks: Diagnosis, Exposure Balance and Correction. Working paper, University of Ljubljana.{p_end}
