# Independent direction screen

This screen concerns Extrusion Net only. It is exploratory, not a registered
confirmation study or a claim of priority. All directions and negative results
will be retained. Selection precedes a larger follow-up study.

## Candidates and decisions

1. **Mechanism identification.** Hold positive loads fixed and observe every raw
   contact. Check whether left/right gates are uniquely recoverable, whether
   observationally equivalent gates disagree under interventions, and whether one
   measured gate intervention resolves the ambiguity. Test both the direct sum
   and the factorized operator. Continue only for an exact structural result,
   beyond numerical optimizer failure.
2. **Saturation as a new operator.** Measure separated-block numerical rank and
   parameter-Jacobian rank before and after Poisson saturation. An invertible
   entrywise transformation cannot fix parameter nonidentifiability. Rank growth
   alone does not establish a new learning principle; compare to established
   nonlinear kernel/feature constructions before selecting this route.
3. **Shared context-dependent gate learning.** Repeat the local-identifiability
   check when a shared logistic gate rule generates all positional gates. This
   guards against incorrectly extending free-gate ambiguity to the actual neural
   controller. One versus multiple contexts tests whether context diversity helps.

Use reproducible float64 numerical checks, finite-dimensional Jacobian ranks,
and seeded heterogeneous loads/gates. Numerical full rank is local evidence,
not a proof of global identifiability. Exclude unused boundary controls. Known
positive loads and known saturation scale are explicit assumptions.

## Follow-up if mechanism ambiguity survives

Derive the exact equivalence family and its admissible range. Test on independent
seeds and lengths. Compare a single intervention at different locations, using
held-out interventions as the outcome. Study conditioning under observation noise;
exact identifiability does not imply stable recovery. Distinguish this inference
problem from the existing network's end-to-end training task.

No candidate receives a high novelty rating merely for passing its experiment.
Priority, practical significance, and competitive validation remain separate gates.
