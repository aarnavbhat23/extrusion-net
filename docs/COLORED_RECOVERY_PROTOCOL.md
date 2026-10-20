# Spaced-probe recovery extension

Theory-led exploratory protocol, written after the four-response strong-barrier
result and before running this extension. Retain that result and its restricted
domain. For any fixed maximum pass probability rho<1, choose K residue classes
so that aliasing from separations K-2, K+2, 2K-2, ... has derivative norm below
one after local two-loading inversion. Prove a projected fixed-point recovery
guarantee and an a posteriori error bound. This is a source-controlled extension
of structured matrix probing, not a new general optimization method.

Use two alternating loading patterns and K one-hot residue-class value vectors,
so the total response count is 2K. When K reaches the number of physical positions,
each color isolates a column and aliasing vanishes. Count vector responses,
observed scalars and input energy; do not call this a four-response result.

Test rho=.2/.5/.8/.95 at 32/128/512 loads, five mechanisms per setting, noiseless
and Gaussian SD .001 outputs, a fixed-point tolerance 1e-10 and at most 100
iterations. Choose K from the analytical bound to give contraction <=.5 before
seeing any mechanism. Report actual gate error, residual certificate, iteration
count, finite/infinite-domain contraction bounds and fit failure. Check the
matrix-free implementation against the production scan and dense response.
Use a new held-out load pattern and one right-gate cut to test response transfer.

Seed 20261111 is an arbitrary reproducibility identifier. A uniform gate bound,
known loading and source-invariant transport remain required. The aim is a global
recovery theorem with a transparent response-budget tradeoff, not a proof of
generality from task labels or physical manipulability of genomic sources.

After the random-mechanism screen: add a near-transparent deterministic stress
case, all gates .999*rho, at rho .8/.95 and lengths 128/512, with zero/.001
noise. This checks a regime where long products are much larger than under
independent uniform gates. Preserve the original random results and response
budgets.

Then run an explicit source-invariance failure check, ten mechanisms each at
rho .2/.5/.8, length 128, sigma=.001. Two calibration loads use only .1/.9.
Construct a gate shift gamma*mean(1-((q-.5)/.4)^2), gamma=0/.005/.02/.05,
clipped into the gate box. It vanishes on both calibration patterns but acts
on a held-out uniform-.5 loading. This is a deliberately adversarial alternative,
not a model of biological perturbations. Use a 97.5% calibration gate-error box
and a 97.5% held-out noise bound to test the new positive-value response. Their
union bound gives at least 95% non-rejection under the invariant null. Report
false rejection, detection and intervals, including undetected shifts. No finite
calibration set can establish invariance under every untested loading rule.
