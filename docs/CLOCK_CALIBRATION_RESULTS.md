# Finite-data clock calibration: failed continuation gate

Protocol: [CLOCK_CALIBRATION_PROTOCOL.md](CLOCK_CALIBRATION_PROTOCOL.md),
committed as `60cf44c` before outcome evaluation. All three mechanisms, two
initiated-episode budgets, four repeats and three strategies are retained:
72 fits, with three fixed optimization starts each. No cases were removed.

## What the experiment establishes

At a fixed number of initiated episodes, retaining absences gives much more
accurate intervention predictions in these matched-family simulations than
discarding them and collecting a second normalized histogram. Widely separated
clocks improve on nearly identical clocks, but fail the frozen precision gate.

The target is the full categorical snapshot at gamma=.3 after doubling removal.
Errors below are medians over all 12 mechanism/repeat cases at each budget.
Absence errors are percentage points; total variation is on the 0–1 scale.

| Initiated episodes | Observation strategy | Absence error (points) | Total variation |
|---:|---|---:|---:|
| 5,000 | Two clocks, .1/.3 | 5.515 | .05515 |
| 5,000 | Close clocks, .28/.3 | 26.511 | .26511 |
| 5,000 | One clock, absences retained | .530 | .00685 |
| 20,000 | Two clocks, .1/.3 | 2.018 | .02018 |
| 20,000 | Close clocks, .28/.3 | 8.357 | .08357 |
| 20,000 | One clock, absences retained | .118 | .00298 |

The primary two-clock medians by mechanism at 20,000 episodes are **2.712,
3.298 and 1.350 points**. Only 1/3 meets the <=2-point tolerance; the protocol
required at least 2/3. The continuation gate **fails**. The small number of
repeats supports a development decision, not a population-level superiority test.

All 216 optimization starts report success. That is a solver termination status,
not proof of global optimality. The selected result has the greatest saved
likelihood among the three starts. Fitting parameters never receive the truth.

All **48** simultaneous 95% two-clock speed confidence intervals have an
unbounded upper endpoint, including all 24 at 20,000 episodes. All contain the
generating speed, which is weak evidence when they are so broad. These intervals
use only the two early-state ratios. Their failure does not prove that every
full-likelihood confidence procedure would fail, and they do not certify the
point predictions in the table.

## Resource and observation accounting

Two-clock likelihoods condition on presence and discard absences and initiation
counts. The full-retention control uses those counts. Thus this compares useful
observation contracts, not competing algorithms given identical information.
The [theory note](NORMALIZED_SNAPSHOT_THEORY.md) states precisely when the
normalization ambiguity applies and proves one-clock recovery when the absolute
denominator is available.

At 20,000 initiated episodes the expected sum of sampled observation durations
is 133,333 time units for .1/.3, 69,048 for .28/.3, and 66,667 for full retention.
These are sums of nominal exponential delays, not measured microscope time,
parallel wall time or actual assay cost. Two separated clocks have twice the
nominal observation duration of the full-retention control at this count budget.
Individual timestamp likelihoods have not been tested.

## Verification and reproduction

Data are generated with a 324-state generator/resolvent; fitting uses an acyclic
flow implementation. Both assume the same capture family. A third sparse-matrix
implementation in `verify_clock_calibration.py` imports neither implementation.
It independently regenerates all saved random counts, checks all 216 start
likelihoods, selected fits, intervention distributions, resource ledgers, interval
endpoints and the failed gate. It does not refit or establish model validity.

Saved observations and fits: [clock_calibration.json](../results/clock_calibration.json).
Independent replay: [clock_calibration_verification.json](../results/clock_calibration_verification.json).

```sh
python -m unittest discover -s tests -v
python verify_clock_calibration.py
# Optional, slower: regenerate all optimization fits.
python study_clock_calibration.py
```

No biological cohort was fitted in this round. Public distance trajectories
still lack the validated known-loading/anchor-state bridge required here.
Computational independence is not biological replication.

## Research decision

Stop promoting two-clock acquisition as the preferred method. Preserve the
exact observation-contract result: normalization can erase intervention-relevant
information; measuring the denominator can recover it. A short theory working
draft is justified as a vehicle for expert scrutiny. Submission as an 8–9-level
paper is not justified yet: priority, timestamp controls, uncertainty in loading
and informative dropout remain unresolved. More fitted variants of this same
family are not the next priority.
