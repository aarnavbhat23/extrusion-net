# Finite-data calibration screen

Freeze before evaluation on 2026-10-28. Extend the capture family to unknown
positive arm speeds. Parameters are (vL,vR,mu,pL,pR,h); known loading and barriers
remain as before. The proposed normalized-resolvent ambiguity and two-clock
inverse are algebraic statements, not a validated biological observation model.

Use three fixed new mechanisms:
(1.15,.85,.27,.35,.65,.45), (.75,1.25,.5,.6,.3,.9),
(1.3,.7,.7,.45,.55,.6). Use total initiated-episode budgets 5,000 and 20,000,
four independent repeats each, seed 20261029. Retain all 24 cases.

Compare three observation strategies at equal initiated-episode counts:
(1) two normalized present-state histograms at gamma=.1 and .3, half per assay;
(2) two close clocks .28 and .3, half per assay, a conditioning control;
(3) one full categorical snapshot at gamma=.3, retaining the absent category.
Retained-state likelihoods do not use the absent category or total initiation
count as data. The resource ledger does report initiated and retained counts and
expected observation-clock duration. Equal episode counts are not equal time or
instrument cost; no economic or timing advantage is implied.

Use full-data multinomial likelihood with three fixed starts:
(1,1,.3,.4,.6,.6), (.6,1.4,.6,.7,.3,1), (1.4,.6,.15,.2,.8,.3).
Domain: speeds [.2,3], mu [.03,2], capture probabilities [.001,.999], h [.1,1.5].
L-BFGS-B, maxiter=300, normalized negative log likelihood, ftol=1e-11.
Save all starts and statuses. The generating parameters are not fit inputs.
Predict the full snapshot distribution at gamma=.3 after doubling removal;
report total variation and absolute absence error. Point predictions are not
confidence certificates. Do not infer superiority from four repeats per mechanism.

For the two-clock strategies, construct simultaneous 95% Clopper-Pearson odds
bounds for q[0,0]/q[1,0], then invert their difference to bound vL. Report unbounded
intervals explicitly. They provide a finite-sample conditioning diagnostic, not
a prediction interval for the intervention. Check their coverage on all cases.

Generate counts using a separately assembled finite-state generator/resolvent,
and fit with an acyclic forward flow. They describe the same assumed family;
this is independent computational verification, not experimental replication.
Two clocks proceed only if median absence error is <=.02 in at least two of the
three mechanisms at N=20,000. If the full-retention control matches or wins,
favor it whenever that measurement is accessible. No policy selection after
these development outcomes will be described as fresh confirmation.

Retained individual timestamps remain a required future control before claiming
measurement optimality. They require a different likelihood and are not silently
represented by the coarse two-clock histograms. Public live-cell distance data
do not currently establish the controlled loading/anchor observation contract;
the biological validation gate remains open.
