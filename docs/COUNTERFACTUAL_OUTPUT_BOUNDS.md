# Sharp bounds for an intervened routing mixer's outputs

A fixed single-gate intervention has more structure than entrywise contact
monotonicity alone suggests. The entire unsaturated counterfactual contact matrix
traces a line segment across the free-gate equivalence family. Consequently signed
linear message bounds are sharp at its endpoints. The production mixer's log-mass
feature needs one additional possible stationary point per output.

This strengthens the earlier result restricted to nonnegative weighted sums. It
is a statement about an exact observed contact map and a specific computational
intervention, not statistical coverage for a noisy map or identifiability of the
shared neural controller.

## Common fractional-linear parameter

Use the notation of `IDENTIFIABILITY.md`, with initial right pass u=r[0]. Write

    P_0(u)=1, P_1(u)=u,
    P_(k+1)(u)=A_k P_k(u)-B_k P_(k-1)(u).

Each P_k is affine in u, and r_k=P_(k+1)/P_k. For a fixed intervention that
multiplies r_k by rho in [0,1], the changed reduced contacts satisfy

    C_rho[a,b](u) = C[a,b] - (1-rho) C[a,k] P_b(u)/P_k(u),  a<=k<b.

Other entries remain fixed. All changed entries share the same denominator P_k.
This remains true after restoring the symmetric contact matrix D. If u runs over
the feasible interval [l,h], then

    D_rho(u) = (1-theta(u)) D_rho(l) + theta(u) D_rho(h),
    theta(u) = ((u-l)/(h-l)) P_k(h)/P_k(u).

With positive interior products, P_k stays positive on the relevant closed
interval. Theta is continuous, increasing, and maps the interval onto [0,1]. The
endpoints describe the closure; some gate probabilities at them may be zero or
one. The resulting sharp extrema may therefore be limits of strictly interior
factorizations. Known positive contact scaling beta preserves the segment.

Every signed linear query of D, including every value-channel message D V, is
affine in theta. Its sharp minimum and maximum are attained at the endpoints
(or approached at the closure). A positive row-normalized message is fractional
linear and also lies on its row's endpoint segment. The current production mixer
uses unnormalized messages plus log mass, so that separate normalized variant is
not substituted for the implemented model.

## The actual one-head mixer

For fixed values V and affine output weights, one output coordinate in the
unsaturated production formula is

    f(theta) = a+b theta+c log(d+e theta),  d+e theta>0.

The linear message produces a+b theta. The mass feature log(1+sum_j D_ij)
produces the logarithm. Its derivative is

    f'(theta)=b+c e/(d+e theta).

If b and e are nonzero, the only possible interior stationary point is

    theta* = -c/b-d/e.

Evaluate theta=0, theta=1, and theta* when it lies inside (0,1). These candidates
give sharp componentwise bounds. Zero b/e cases are constant or monotone and need
only endpoints. Classification margins formed by an affine readout have the same
form and can be bounded using the corresponding difference of output weights.
This does not certify a subsequent arbitrary nonlinear classifier.

`counterfactual_bounds.py` implements this result. For the message-only part,
two endpoint scan evaluations suffice. Given those endpoint messages and masses,
the log-mass correction uses scalar operations per output. The demonstration
constructs dense endpoints for validation; it does not claim a new fused kernel.
The existing linear scan can evaluate the endpoints without materializing dense
contact matrices when the endpoint gates and values are already available.

## Validation and a necessary interior extremum

An exploratory screen used 60 mechanisms, matrix sizes 8/16/32, random signed
values and affine output weights, and fixed gate multipliers 0, .3 or .7. It
checked 101 family points per mechanism, then reconstructed physical gate values
at every claimed extremum and evaluated the independent original contact sum.

- Maximum matrix-segment discrepancy: 1.47e-13.
- Maximum output enclosure discrepancy: 7.11e-15.
- Maximum extremum-attainment discrepancy: 7.11e-15.
- Endpoint-only output bounds missed 24 of 7,080 output coordinates.

One retained witness has endpoint output range [6.29155,6.34611], but reaches
6.42563 at theta=0.46181. Thus blindly propagating only endpoint *network outputs*
is incorrect when the mass feature is included. The three-candidate rule captures
that extremum. Unit tests also compare against the production scan and independent
scalar optimization, including signed values, constant intervals and rejection
of an inconsistent noisy baseline map.

These are numerical checks of a real-arithmetic derivation. The implementation
is not outward-rounded interval arithmetic, and near-degenerate or underflowing
families require separate numerical treatment. A noisy map projected onto a
feasible one does not turn these bounds into a statistical confidence interval.

## Research scope and priority

The result applies to one fixed directional-gate multiplier and one unsaturated
head with fixed values and affine output. Multiple simultaneous interventions,
Poisson saturation, unknown values, arbitrary nonlinear downstream layers and
shared-controller constraints are not covered by this implementation. Applying
free-gate bounds to a more constrained controller is potentially conservative;
sharpness is established for the stated free-gate family, not every neural model.

Fractional-linear optimization and partial-identification bounds are established;
for example [counterfactual policies using instruments](https://www.jonathandroth.com/assets/files/EvaluatingCounterfactualsWithIV.pdf)
uses fractional-linear programs, and [optimal experiments under partial identification](https://arxiv.org/abs/2605.06993)
studies tightening identified sets. The specific common-denominator reduction
for this routing operator is the result investigated here. Priority remains
unproven. This supplies an exact restricted inference primitive; it is not yet
evidence of a better general learning architecture.
