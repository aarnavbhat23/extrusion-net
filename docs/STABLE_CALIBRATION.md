# Identifiability is not stable recovery

The exact free-gate theorem says that one scalar calibration can remove a
one-dimensional ambiguity. Noisy-map experiments show that this is not a guarantee
of accurate recovery. The following local calculation explains a distinct issue:
weakly observed directions can remain after the exact null direction is removed.

## Model and scope

Let there be p directional gate pairs and p+1 loads, with q=1 and homogeneous
interior passes l=r=a, 0<a<1. Consider the local submodel

    l_i = a + z_i,    r_i = a - z_i.

Thus loads and the sums l+r are known. This is an easier oracle problem than
estimating every load and gate. Observe every upper contact with independent
Gaussian noise of variance sigma^2. A calibration at site s measures the local
right-cut response C_cut[s,s+1]=a+z_s with independent noise of the same variance.
All results below concern Fisher information at z=0 and locally unbiased
estimation. They do not bound every biased estimator, impose a prior on a neural
controller, or establish global finite-noise minimax risk.

## Contact Jacobian and spectral bounds

For a contact spanning d gates, starting at i, differentiation gives

    delta C[i,i+d] = a^(d-1) sum_{r=0}^{d-1} (d-1-2r) z[i+r]
                  = a^(d-1) sum_{r=1}^{d-1} r(d-r)
                    (z[i+r-1] - z[i+r]).

The d=1 derivative is zero. At d=2 it is a(z_i-z_{i+1}). Let J contain all
contact derivatives, F=J'J, and L be the unweighted path Laplacian on p gates.
The d=2 observations alone imply F >= a^2 L in positive-semidefinite order.

For each d use weighted Cauchy-Schwarz, then sum over all starting positions.
Each difference is counted with total weight at most
W_d=sum_{r=1}^{d-1} r(d-r)=d(d^2-1)/6. Therefore

    a^2 L <= F <= c(a) L,
    c(a) = sum_{d=2}^infinity a^(2d-2) W_d^2
         = r(1+r)(1+8r+r^2)/(1-r)^7,  r=a^2.

The finite sum is bounded by this convergent infinite sum. The closed form follows
from the power-series sums of d^2, d^4 and d^6. The uniform z direction is the exact
null direction; slowly varying directions have small information as well.

## One boundary calibration

Adding the boundary observation gives information
(F+e_0 e_0')/sigma^2. Inverting the spectral bounds and using the anchored path
inverse, whose diagonal is 1+i/c for cL+e_0e_0', gives

    sigma^2 [1+(p-1)/(2c(a))]
      <= trace(CRLB)/p
      <= sigma^2 [1+(p-1)/(2a^2)].

For every fixed 0<a<1, this average local variance bound grows linearly with
chain length. A scalar that is sufficient for exact identifiability is therefore
not sufficient for a length-independent local recovery guarantee.

## Any k calibration sites

Give an oracle the exact values at the k sites, making recovery easier still.
The remaining p-k sites form at most k+1 unanchored segments. A segment of length
ell between two exact anchors has grounded-path inverse trace
ell(ell+2)/(6c); an end segment has trace ell(ell+1)/(2c). Each is at least
ell^2/(6c). Cauchy-Schwarz over segment lengths yields

    trace(CRLB)/p >= sigma^2 (p-k)^2 / [6 c(a) p (k+1)].

This applies to any placement of k>=1 measured calibration sites under the
stated oracle submodel. It is not an assertion that all placements attain it.
Conversely, a regularly spaced set of noisy anchors and the d=2 differences
suffice for local variance at most
sigma^2 [1+ceil(p/k)/a^2]: reconstruct each site from its nearest anchor along
the path. Thus the order p/k is achievable in the local linear experiment.
Holding normalized local variance bounded as p grows requires a calibration
budget proportional to p, for fixed a, under these assumptions.

This is a special-case stability result, not a universal impossibility result
for Extrusion Net. Shared controllers and informative priors can reduce the
number of degrees of freedom. At high a the constants are loose; numerical
bounds must not be advertised as observed nonlinear errors.

## Relationship to existing research

Laplacian information bounds and anchor placement are established in
[rotation synchronization](https://arxiv.org/abs/1211.1621) and
[localization](https://www.ocf.berkeley.edu/~chechang/docs/Eurasip06.pdf).
The contribution investigated here is the explicit reduction and calibration
scaling for this extrusion operator. Priority for this specialization has not
been established. The local design rules are classical experimental design.

`study_calibration_budget.py` checks the spectral inequalities at 18 size/pass
settings, compares calibration placements and budgets, and simulates the
linearized Gaussian experiment. `tests/test_calibration_budget.py` checks the
analytic derivative against the independent nonlinear contact sum, the closed
series, the spectral inequalities and the anchored inverse. None of those tests
substitute for a useful learned routing model on an independent task.
