# Common observable and intervention screen — frozen before outcomes

Use an exponentially randomized observation delay T~Exp(gamma) after controlled
loading of one motor at site 8 on a 17-site lattice. The observation is its
anchor pair if it survives until T, or an absent category if removed earlier.
This is an experimental-design toy observable, not Hi-C or polymer contacts.
Each independent episode gives one categorical observation, avoiding unbounded
residence-time measurements. Use unmodified lefsim 0.1.0 and preserve provenance.

Establish the exact resolvent bridge: for a position-Markov process with left
rates L_i, right rates R_j and constant removal mu, snapshot mass obeys

    C_ii = gamma*q_i/(L_i+R_i+mu+gamma)
    C_ij = (L_{i+1} C_{i+1,j} + R_{j-1} C_{i,j-1}) /
           (L_i+R_j+mu+gamma), i<j.

With the corresponding index and constant mappings this must agree with the
unmodified authors' netdlem2 forward, transform=False, constant=1. Source:
chikinalab/dLEMpytorch commit e55f9226192c34c34990518713a15ec42dff884a.
Download the source outside tracked code, retain its hash and attribution, and
do not redistribute it under this repository's license. Validate the bridge to
1e-10 on independent small-rate examples before comparing predictions.

Develop a finite-state capture model that retains whether either arm has been
captured. Capture occurs on arrival at sites 5/11; captured arms cannot step.
The removal rate is mu*h if either is captured, mu otherwise, including bounds.
Compute categorical snapshot probabilities by an acyclic flow/resolvent solve.
This is classical Markov modeling of an existing engine mechanism, not a new
biophysical process or architecture. An independent engine checks its semantics.

Point-comparison cases: mu=.15/.35/.7, (pL,pR)=(.2,.8)/(.5,.5)/(.8,.2),
h=.25/1: 18. Add three off-grid cases (.25,.35,.65,.6), (.55,.7,.4,.4),
(.9,.4,.7,.8). At gamma=.3 acquire 20,000 baseline observations and 20,000
independent observations after multiplying all removal rates by two. Compare
the continuous four-parameter capture fit, a finite-bank point fit, and the
published dLEM recurrence with all 16 effective step rates plus removal fitted.
All use the same baseline categories and known source. Fit by squared Hellinger
distance with three fixed initializations; retain convergence and training errors.
Validate the fitted NumPy recurrence against the unchanged authors' forward;
comparative outputs use that forward. Report the extra parameters, numerical
optimization limitations, time, all cases and target RMSE. Continue the hidden-
state prediction claim only if median target RMSE improves at least twofold
over the strongest dLEM fit, without a median baseline-fit error advantage
exceeding twofold. Otherwise retain the negative result and narrow the claim.

Separate finite-bank acquisition diagnostic: mu=.15/.35/.7, pL,pR=.2/.5/.8,
h=.25/1, giving 54 supplied candidate mechanisms. Use the 18 on-grid mechanisms,
three seeds each. Actions choose gamma=.1/.3/1/3; each adds 1,000 observations,
up to eight batches. Aggregate counts per action. Clopper-Pearson intervals use
alpha=.05/(4*8*82) per category and sampling horizon. Intersect candidate sets.
Under membership in the finite bank this gives simultaneous enclosure across
adaptive choices; it is not a continuous-parameter confidence set. Empty sets
reject the bank. Predict the doubled-removal gamma=.3 distribution; stop when
every target category range is <=.02, otherwise abstain after the budget.

Compare round-robin actions, global Hellinger disagreement per quoted cost,
and target-relevant pair disagreement per quoted cost. The latter prioritizes
candidate pairs whose target distributions differ by more than .02. Cost per
episode is the hypothetical upper bound 1+1/gamma, not a lab timing estimate.
Require >=20% mean cost reduction over the stronger control with no decrease
in precision attainment; all outcomes and abstentions remain visible. Test the
three off-grid mechanisms separately, with no membership/coverage guarantee.
Do not tune action scores after outcomes. Use disjoint seed blocks for every
case, method, action stream and independent target.
