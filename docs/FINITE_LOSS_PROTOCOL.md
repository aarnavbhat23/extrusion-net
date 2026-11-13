# Finite competing-loss prediction screen

Freeze 2026-10-29 before outcome evaluation. Parameters (vL,vR,mu,pL,pR,h,delta).
Use three fresh mechanisms (1.05,.95,.32,.4,.7,.55,.12),
(.85,1.15,.55,.65,.25,.8,.2), (1.25,.75,.65,.5,.6,.4,.08).
Initiated counts 6000 and 24000, three repeats, seed 20261030.

Compare (1) normalized snapshots at (f,gamma)=(1,.1),(1,.3),(2,.3);
(2) absolute snapshots including absence at (1,.3),(2,.3);
(3) absolute snapshots at known scheduled times .5,1,2,4 for each f=1,2.
Allocate the same total count equally across each design's conditions. The third
control retains exact scheduled timestamps, not individual random exponential
timestamps. All observation contracts and nominal durations are reported;
equal counts are not equal acquisition costs or equivalent information.

Split each condition's initiated count in half before drawing independent
multinomial calibration counts. Fit only the first half with three fixed starts:
(1,1,.3,.4,.6,.6,.15), (.6,1.4,.6,.7,.3,1,.3),
(1.4,.6,.15,.2,.8,.3,.05). Domain: speeds [.2,2.5], mu [.03,1.5],
captures [.001,.999], h [.1,1.5], delta [.001,1]. L-BFGS-B max250,
normalized NLL, ftol=1e-10. Save all statuses and counts; truth is never a fit input.

The second half defines the conventional split-likelihood confidence set
L_holdout(theta)>=.05*L_holdout(theta_hat_train). Check generating-parameter
membership in matched simulations. This is not a predictive interval or a
proof of useful set width. Do not turn finite optimizer searches into certified
outer intervals. Small repeat counts cannot validate 95% coverage empirically.

Predict the full snapshot at f=3,gamma=.3. Report absence error and total
variation. Primary point-accuracy gate: median absence error <=.02 in >=2/3
mechanisms at N=24000; this is a computational screen, not a clinical tolerance.
Report each method separately and keep wins by simpler controls.

Stress cases: N=24000, first two repeats, all three mechanisms/designs;
(a) actual calibration f=2 means 1.8 but fitting assumes 2; target f=3 is true 3;
(b) extra loss .15 once either arm reaches >=4, in calibration and target,
while fitting assumes state-independent delta. Keep all 90 fits. These stresses
invalidate the matched-model confidence theorem; membership is only diagnostic.

Data use an independently assembled 324-state generator (resolvent for random
clocks, matrix exponential for scheduled snapshots). Fitting uses the existing
flow recurrence and a 100-state reachable-state uniformization. Verify both
backends before evaluation. Uniformization uses fixed rate20 and a Poisson tail
below 1e-13; certify it covers every fitted exit rate. Solver convergence is not
global optimality. No biological validation or novelty upgrade follows alone.
