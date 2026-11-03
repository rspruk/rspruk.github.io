{smcl}
{* *! version 0.9.6}{...}
{title:Title}

{phang}
{bf:splitsc} {hline 2} Time-split honest inference for synthetic controls with imperfect fit

{title:Syntax}

{p 8 16 2}
{cmd:splitsc} {depvar} {ifin}{cmd:,}
{opt u:nit(varname)} {opt t:ime(varname)} {opt tr:eated(string)} {opt trp:eriod(#)} {opt train:end(#)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt u:nit(varname)}}unit identifier, string or numeric{p_end}
{synopt:{opt t:ime(varname)}}numeric time variable{p_end}
{synopt:{opt tr:eated(string)}}value of {it:unit} for the treated unit, a number if {it:unit} is numeric{p_end}
{synopt:{opt trp:eriod(#)}}first post-treatment period{p_end}
{synopt:{opt train:end(#)}}last period of the training block. Periods after it and before {opt trperiod()} form the validation block{p_end}

{syntab:Options}
{synopt:{opt noint:ercept}}simplex weights without an intercept. The default is de-meaned synthetic control (simplex weights with an intercept){p_end}
{synopt:{opt exc:lude(numlist)}}periods dropped from every block, for example years with possible anticipation{p_end}
{synopt:{opt c:urvature(#)}}also report the interval under the curvature bound {it:C} = {it:#}{p_end}
{synopt:{opt l:evel(#)}}confidence level, default {cmd:c(level)}{p_end}
{synopt:{opt phi:bar(numlist)}}grid for the breakdown persistence, values in [0,1), default 0.1(0.1)0.9{p_end}
{synopt:{opt rep:s(#)}}simulation draws for the AR(1) and Theorem 3 critical values, default 20000. Use 200000 or more for reported results{p_end}
{synopt:{opt seed(string)}}random-number seed{p_end}
{synopt:{opt exch}}also report the interval with the placebo-calibrated drift bound (Theorem 3). Cannot be combined with {opt noplacebo}{p_end}
{synopt:{opt nopl:acebo}}skip the placebo fits (no {it:C_don}, no placebo count){p_end}
{synopt:{opt st:ep(#)}}grid step for the Peano-kernel integral, default 0.25 time units{p_end}
{synopt:{opt gr:aph(list)}}draw figures. {it:list} is one or more of {cmd:paths}, {cmd:gap}, {cmd:sensitivity}, {cmd:placebo} and {cmd:persistence}, or {cmd:all}{p_end}
{synopt:{opt gn:ame(prefix)}}prefix for the names of the graphs, default {cmd:splitsc_}{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:splitsc} fits synthetic-control weights on a training block of pre-treatment periods only, and uses
a later validation block to estimate the average post-treatment effect with the trend-balancing anchor.
The anchor extrapolates a line fitted to the validation-block gap to the mean post-treatment period, so the
linear-class interval is the classical prediction interval for the mean of the post-treatment gaps. Because the
weights do not use the validation or post-treatment outcomes, the interval has finite-sample coverage under
Gaussian, serially independent noise, whatever the fitting rule.

{pstd}
The curvature class bounds the second derivative of the counterfactual discrepancy by {it:C} and leaves its
slope free. The bias bound is {it:C} times the curvature exposure {it:kappa_C} of the anchor, computed by a
Peano-kernel integral. The breakdown curvature {it:C*} is the smallest {it:C} at which the interval includes
zero. Placebo fits of every donor against the other donors (treated unit excluded) give the placebo
curvatures, their maximum {it:C_don} and the number of placebos whose curvature exceeds {it:C*}. These are a
sensitivity analysis, since they use donor outcomes after the training block.

{pstd}
The breakdown persistence is the smallest value on the {opt phibar()} grid at which the interval, with the
Student quantile replaced by a numerical critical value for stationary AR(1) disturbances, includes zero. For
each value of phibar the critical value is the largest simulated quantile over eleven equally spaced
coefficients from 0 to phibar. The calibration assumes that the dependence does not reach back into the
training block.

{pstd}
With {opt exch}, the drift bound is calibrated from the placebo errors, inflated by the treated unit's own
noise, and the anchor is chosen from a one-parameter family with a simultaneous critical value. Coverage is at
least 1 - alpha - 1/(J+1) if the treated unit's placebo error is exchangeable with the donors'.

{title:Figures}

{pstd}
{opt graph()} draws up to five figures, each kept in memory under the name {it:prefix}{it:figure}, for example
{cmd:splitsc_gap}. Save or export them with {helpb graph export}.

{phang2}{cmd:paths}: the treated unit and its synthetic control, with the validation and post-treatment blocks
shaded. The synthetic control includes the intercept when the weights are de-meaned.{p_end}
{phang2}{cmd:gap}: the gap, the validation-block trend extrapolated over the post-treatment block, the
post-treatment mean gap (diamond), and the confidence band for the counterfactual mean, centred on the trend at
the mean post-treatment period. The estimate is the vertical distance from the diamond to the centre of the band,
and the interval excludes zero when the diamond lies outside the band. With placebos, a wider band shows the
interval at {it:C_don}.{p_end}
{phang2}{cmd:sensitivity}: the interval as a function of the curvature bound {it:C}, with {it:C*}, {it:C_don},
the bound given in {opt curvature()} and the placebo curvatures marked.{p_end}
{phang2}{cmd:placebo}: the placebo curvature of each donor against {it:C*}. Not available with
{opt noplacebo}.{p_end}
{phang2}{cmd:persistence}: the interval as a function of phibar over the {opt phibar()} grid, with phibar*
marked. This figure needs the critical value at every grid point, so it adds computing time. The extra draws
are made after all other results, which are therefore the same with and without the figure.{p_end}

{title:Requirements}

{pstd}
A balanced panel in long form, with one nonmissing outcome per unit and period in the estimation sample. At
least three validation periods and two donors. The curvature-optimal anchor of the paper, which needs a
second-order cone programme, is not implemented. The trend-balancing anchor is reported throughout.

{title:Remarks}

{pstd}
When the training block is shorter than the number of donors, the simplex least-squares problem can have
several solutions with the same fit. {cmd:splitsc} uses the Lawson-Hanson algorithm, and a different solver
can return other weights with the same fit. The treated unit's interval is typically unaffected, but a
placebo's curvature, and so the count of placebos above {it:C*}, can change by one.

{pstd}
Import outcome data in double precision ({cmd:import delimited ..., asdouble}) to reproduce published values
to more than five digits.

{title:Stored results}

{pstd}{cmd:splitsc} stores the following in {cmd:r()}:{p_end}
{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:r(theta)}}estimate of the average effect{p_end}
{synopt:{cmd:r(se)}}standard error{p_end}
{synopt:{cmd:r(omega)}}residual scale on the validation block{p_end}
{synopt:{cmd:r(nu)}}residual degrees of freedom{p_end}
{synopt:{cmd:r(tcrit)}}Student critical value{p_end}
{synopt:{cmd:r(lb)} {cmd:r(ub)}}linear-class interval{p_end}
{synopt:{cmd:r(C)}}curvature bound from {opt curvature()}{p_end}
{synopt:{cmd:r(lb_C)} {cmd:r(ub_C)}}interval at {opt curvature()}{p_end}
{synopt:{cmd:r(kappa)}}curvature exposure of the anchor{p_end}
{synopt:{cmd:r(Cstar)}}breakdown curvature{p_end}
{synopt:{cmd:r(Cdon)}}largest placebo curvature{p_end}
{synopt:{cmd:r(lb_Cdon)} {cmd:r(ub_Cdon)}}interval at {it:C_don}{p_end}
{synopt:{cmd:r(above)}}placebos with curvature above {it:C*}{p_end}
{synopt:{cmd:r(Pimax)}}largest absolute placebo error{p_end}
{synopt:{cmd:r(phistar)}}breakdown persistence (missing if above the grid){p_end}
{synopt:{cmd:r(phimax)}}largest value of the {opt phibar()} grid{p_end}
{synopt:{cmd:r(ksim)} {cmd:r(Mhat)}}Theorem 3 critical value and drift bound{p_end}
{synopt:{cmd:r(exch_h)}}drift distance of the anchor chosen under Theorem 3{p_end}
{synopt:{cmd:r(exch_theta)}}Theorem 3 estimate{p_end}
{synopt:{cmd:r(exch_lb)} {cmd:r(exch_ub)}}Theorem 3 interval{p_end}
{synopt:{cmd:r(exch_guarantee)}}coverage guarantee 1 - alpha - 1/(J+1){p_end}
{synopt:{cmd:r(intercept)}}intercept of the weights (0 with {opt nointercept}){p_end}
{synopt:{cmd:r(level)}}confidence level{p_end}
{synopt:{cmd:r(J)} {cmd:r(n)} {cmd:r(m)} {cmd:r(n_train)}}donors and block lengths{p_end}
{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:r(cmd)}}{cmd:splitsc}{p_end}
{synopt:{cmd:r(treated)}}treated unit{p_end}
{synopt:{cmd:r(donors)}}donor units, in the order of {cmd:r(weights)}{p_end}
{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:r(weights)}}donor weights{p_end}
{synopt:{cmd:r(anchor)}}weights of the trend-balancing anchor on the validation periods{p_end}
{synopt:{cmd:r(gap)}}time, gap and block indicators{p_end}
{synopt:{cmd:r(placebo)}}placebo errors and placebo curvatures, in the order of {cmd:r(donors)}{p_end}
{synopt:{cmd:r(paths)}}time, treated outcome and synthetic control{p_end}
{synopt:{cmd:r(trend)}}intercept and slope of the validation-block trend of the gap{p_end}
{synopt:{cmd:r(persist)}}phibar, critical value and interval on the grid (with {cmd:graph(persistence)}){p_end}

{title:Example}

{phang2}{cmd:. import delimited using germany.csv, clear asdouble}{p_end}
{phang2}{cmd:. splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exch seed(1)}{p_end}
{phang2}{cmd:. splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) exclude(1990)}{p_end}
{phang2}{cmd:. splitsc y, unit(unit) time(year) treated("West Germany") trperiod(1991) trainend(1975) graph(all)}{p_end}
{phang2}{cmd:. graph export gap.png, name(splitsc_gap) width(2000) replace}{p_end}

{title:Reference}

{pstd}
Honest Inference for Synthetic Controls with Imperfect Fit. Manuscript, 2026.
