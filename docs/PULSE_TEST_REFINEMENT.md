# Post-screen exact-test check — frozen before refined power outcomes

The primary balanced-schedule screen failed its Hoeffding-based power gate.
Do not interpret that as an information-theoretic power limit. Before stopping,
check a stronger conventional exact test on the SAME synthetic mechanisms.

For observed success counts x,y and <=r altered calls per arm, take the supremum
of the equal-arm two-sided Fisher p-value over all compatible true counts.
With equal arm sizes its two-sided p is twice a symmetric hypergeometric tail.
If the compatible count intervals overlap, the supremum is one. Otherwise it
is attained at the nearest interval endpoints. Monotonicity follows from the
hypergeometric success-count coupling when either observed arm count moves
toward the other. Enumerate all count rectangles at small n against this formula.

The supremum p-value is conservative under any allowed corruption, since the
uncorrupted true count pair is among its candidates. This is standard nuisance-
parameter maximization of an exact test, not a new test principle. It is fixed-
time, not valid for selecting a sample size after inspecting these p-values.

Replay seed 2101 and the original generator order to reconstruct the original
Bernoulli count arrays. Verify their ideal Fisher and Hoeffding rejection rates
match every saved row. Apply the corruption adversarially away from the null
for null cases, toward it for alternative cases, as before. Report all schedule
families. The original gate stays failed; this is a development diagnostic of
statistical conservatism, not a fresh confirmation or model-selection run.
