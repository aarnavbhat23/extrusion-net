# Independent event-engine transfer: frozen gate

Use the public lefsim 0.1.0 release, publisher source commit
`db4f5dadd28d16d0c76911202ff9d1b6205bd70e`, without modifying its engine. It is a
new alpha package, not an established benchmark or measured biology. Preserve
its scientific attribution and record provenance. No results from the following
cases have been read before this protocol commit.

An episode starts one motor with both legs at site 8 on a 17-site lattice and
stops immediately before its first Respawn/Evict. Extract outward reach from
the terminal positions. Do not use the new post-respawn position as a terminal
anchor. The earlier bond-origin model must be adapted to a coincident site
origin; this is a documented observable bridge, not direct Hi-C transfer.

Matched engine: step rate 1, all basal/blocked/captured removal rates equal to
mu; convergent CTCFs at sites 5/11 capture with probabilities pL/pR. From this
origin no backside eviction is encountered. Captures stop that arm, the other
continues. Leaving a captured site has pass probability zero, so a probabilistic
capture acts as a probabilistic gate on reaching *beyond* the site. Finite
boundaries clamp reach. Fit mu and both gate strengths using only first-left-step
reach and left/right reach beyond the barrier. Other nonbarrier gates are known
one; the gate locations are given. This is a deliberately restricted transfer.

mu=.15/.5/1, (pL,pR)=(.2,.8)/(.8,.2), seeds 0/1/2: 18 cases. Each uses 10,000
independent calibration episodes and 10,000 new episodes for each uncut/left
capture=1/both captures=1 target. Point predictions compare asynchronous,
synchronized and independent-arm kernels with the same fitted marginal survival
and gates. All 81 joint-reach events are targets. Account for binomial calibration
uncertainty using simultaneous Clopper-Pearson intervals on the three fitted
probabilities, and propagate monotone parameter boxes. Independently bound target
Monte Carlo noise by a Hoeffding union bound. Alpha .025 calibration + .025 test
per case across all three targets. No useful-width guarantee is assumed.

Negative controls: three mu values, pL=pR=.8, one fixed seed; either captured
removal slowed to mu/4, or right CTCF flipped to evict from inside. Fit the same
restricted model and test it; retain missed mismatches. Full episode samples are
not laboratory measurements and repeated trials are not biological replication.

Continue this transfer only if >=17/18 matched cases are not rejected, the mean
point RMSE is at least twofold lower than the independent-arm comparator, and
at least 4/6 mismatch cases reject on new targets. Report all widths and extra
simulation costs. A failed mismatch criterion blocks a broadly valid biological
model claim even if the matched kernel reproduces the simulator.
