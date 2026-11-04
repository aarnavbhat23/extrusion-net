# Competing-loss calibration: structural screen

Freeze on 2026-10-28 before running the numerical screen. Nature motivates two
distinct processes: regulated motor removal and an unaffected observation-loss
process. This is a mathematical model, not an assertion that a particular drug
selectively changes one rate. Generic identifiability analysis is established.

Extend the six-parameter capture family by a common nonnegative loss delta.
At turnover multiplier f and relative speed s, use free removal f*mu+delta,
captured removal f*mu*h+delta, and speeds s*vL,s*vR. Gamma is known. Capture
probabilities, loading and geometry are unchanged across conditions.

Use all three mechanisms from the preceding clock study, each with delta=.15.
Evaluate six declared designs: relative speed (.5,1,2) at f=1,gamma=.3;
turnover (.5,1,2) at s=1,gamma=.3 with delta known zero; the same turnover with
unknown delta; two clocks (.1,.3) at f=1 with unknown delta; three conditions
(f,gamma)=(1,.1),(1,.3),(2,.3) with unknown delta; absolute-retention observations
at (1,.3),(2,.3) with unknown delta. Normalize only the first five designs.

Verify matrix-resolvent probabilities against the separate flow implementation.
Construct explicit equivalence transformations with scale c=1.6, and report all
condition distances. A failing gauge is not an identification proof. Compute
central-difference sensitivity singular values at relative steps 1e-4 and 5e-5;
these are local diagnostics, never proofs of global identification.

Derive an explicit inverse for the three-condition and absolute-retention
designs. Check all mechanisms and both step sizes; require maximum reconstruction
error <1e-9. Include one-sided capture and no-capture exceptions in tests.
Report a no-capture h as irrelevant, not identified. Do not fit finite data or
claim acquisition savings in this screen. Continue only to finite-data design
if the inverse checks pass and its measurement assumptions remain explicit.
Unknown intervention strengths and arbitrary dropout/loading remain unresolved.
