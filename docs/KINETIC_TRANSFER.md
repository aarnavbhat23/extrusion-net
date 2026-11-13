# What survives finite motor lifetimes?

The original recovery algorithm transfers to independent arm lifetimes but
fails under a shared removal clock. A local contact-band diagnostic identifies
the clock in the idealized shared-clock family. That family can also have
unbounded off-diagonal rank, so it exceeds the existing fixed-state linear scan.
These are structural results in an independently enumerated simulator, not
validation on measured chromatin or a new task-learning result.
The generator assumes discrete equal arm speeds, permanent gate failure and
geometric removal, without collisions or polymer motion.

## Experiment and retained failures

[Protocol](KINETIC_TRANSFER_PROTOCOL.md), [generator](../study_kinetic_transfer.py),
[all 240 fits and post-screen rank witnesses](../results/kinetic_transfer_study.json).
The protocol was written locally before execution; this is exploratory work,
not a separately committed confirmation. Ten mechanisms (five each at 16/32
loads), four survival probabilities and three output-noise levels were used
for each lifetime model. The original recovery code and analytical color choice
were unchanged. Calibration and held-out Gaussian bounds each used alpha .025.

| Mechanism | Survival | Held-out rejection, SD 0 / .0001 / .001 | Noiseless mean message RMSE |
|---|---:|---:|---:|
| Independent clocks | 1 / .95 / .8 / .5 | 0/10 in every setting | below 6e-16 |
| Shared clock | 1 | 0 / 0 / 0 of 10 | below 6e-16 |
| Shared clock | .95 | 10 / 10 / 0 of 10 | .01820 |
| Shared clock | .8 | 10 / 10 / 10 of 10 | .05103 |
| Shared clock | .5 | 10 / 10 / 10 of 10 | .05727 |

All fits converged. A good marginal-gate fit is insufficient: noiseless
shared-clock gate RMSE versus effective marginal gates was below 6.2e-7,
yet held-out uncut/cut messages disagreed. At SD .001 the conservative validation
interval misses the weaker s=.95 mismatch. Absence of rejection is not evidence
of model correctness. The independent-clock control is algebraically matched
to the original family and therefore supplies no evidence of broader expressivity.

## Conditional identification of a shared clock

Use physical positions 0,...,m and loading q_k at bond k. Gates l_h and r_h
control the left and right reach. A geometric lifetime T with P(T>=t)=s^t
shared by both arms gives

\[
D_{ij}=\sum_{k=i}^{j-1}q_k
 \left(\prod_{h=i}^{k-1}l_h\right)
 \left(\prod_{h=k}^{j-2}r_h\right)
 s^{\max(k-i,j-1-k)}.
\]

Independent lifetimes give exponent j-i-1 instead, exactly absorbing s into
both gate vectors. In the shared-clock case, two known loading patterns with
nonzero local determinants recover a_h=s l_h and b_h=s r_h from the
distance-two contact band, just as in the earlier local inversion. At distance
three, subtracting the static prediction gives

\[
D_{i,i+3}-(q_i b_i b_{i+1}+q_{i+1}a_i b_{i+1}
  +q_{i+2}a_i a_{i+1})
=q_{i+1}a_i b_{i+1}(s^{-1}-1).
\]

Thus any nonzero middle coefficient identifies s, then l=a/s and r=b/s,
**conditional on this model family**. Across locations, the inferred ratio
must agree; disagreement falsifies the family. On clean maps the recovered
s agrees to numerical precision (ratio spread below 4e-14). These extra
contact bands are not part of the original vector-response budget. No noisy
clock-recovery guarantee or biological parameter estimate is claimed. Under
independent clocks this diagnostic reports s=1; it cannot identify each arm's
physical lifetime separately from its pass gates.

## Exact rank obstruction

Place one source at a central bond, with r positions on either side and all
pass gates one. Index each side by distance a,b=0,...,r-1 from that bond.
The cross-block is H_ab=s^max(a,b). For 0<s<1,

\[
\det H=s^{r(r-1)/2}(1-s)^{r-1}>0.
\]

Proof: H=U diag(d) U^T, with U_at=1[a<=t],
d_t=s^t(1-s) for t<r-1 and d_(r-1)=s^(r-1).
U is triangular with unit diagonal, giving the determinant above. Hence the
off-diagonal block has rank r for arbitrary r. Positive gates merely introduce
invertible row/column scaling. By continuity, sufficiently small positive
loads at other bonds preserve full rank, so zero loading is not essential.

A fixed-dimension **linear** sequential state can transmit at most its state
dimension across a sequence cut. It therefore cannot exactly implement every
member of this shared-clock family. This does not rule out other linear-time
algorithms, nonlinear computational models, or accurate finite-rank
approximations. Numerical witnesses at r=4/8/16/32 and s=.8 are supplementary;
the determinant proof, not numerical rank thresholding, establishes the result.
The max kernel and rank argument use elementary established matrix facts;
their application here is not a claim of a new general matrix theorem.

## Decision and next gates

Continue the mechanism-identification direction. Do not present the static
operator as a general kinetic model. The next bounded experiment should freeze
clock-aware inversion on noisy local bands, assess clock/gate uncertainty and
predict new loads plus two-sided cuts generated by a separate kinetic engine.
Compare directly with the original static model and a flexible fitted contact
baseline; count all observed scalars. Stop the extension if clock estimates
are unstable at useful noise levels or if mismatch cannot be detected.

Only after that gate should approximation rank versus error and end-to-end task
value be studied. A measured-data transfer requires independently known source
and intervention information; fitting an observational contact map alone does
not test the current controlled-probe theorem.

Finite motor residence time/processivity already appears in
[Gassler et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5730859/).
The contribution being investigated is recoverability, falsifiability and
computational consequences under explicit probe access, not finite lifetime
as a new biological mechanism.
