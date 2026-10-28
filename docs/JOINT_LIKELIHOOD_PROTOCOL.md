# Joint likelihood audit, frozen 2026-10-28

Question: do the earlier rectangular-calibration ambiguity witnesses survive a
joint likelihood confidence set? This is a development reanalysis of the same
six already evaluated datasets, not a new confirmation cohort. No biological
claim or new inference method is proposed.

Use all six rows of `results/continuous_occupancy_study.json`, in their saved
order. Split each 10,000-count assay into two equal parts by multivariate
hypergeometric sampling with one generator seeded 20261028. Save every split.
Fit the four-parameter capture family jointly to the two training assays by
multinomial maximum likelihood, bounded by the existing continuous domain.
Use fixed starts (.15,.2,.2,.25), (.35,.5,.5,.6), (.7,.8,.8,1.), L-BFGS-B,
at most 300 iterations per start. Generating parameters are not fitting inputs.

Let q be the fitted training distribution, L the product likelihood on the two
held-out assays, and alpha=.05. Retain theta when

    log L(theta) >= log L(q) + log(alpha).

Conditional on training, L(q)/L(theta_true) has expectation at most one.
Markov's inequality gives at least 95% coverage per dataset under the specified
family and independent sampling. This is the split likelihood construction of
[Wasserman, Ramdas and Balakrishnan](https://arxiv.org/abs/1912.11436), not a new
guarantee. Reanalysis after seeing data is exploratory; nominal coverage is a
property of a prospectively fixed procedure, not a retrospective validation.

Report held-out likelihood slack for all twelve old rectangular witnesses and
all six generating parameters. Search for new low/high intervention-absence
witnesses using SLSQP, the fitted training parameter and the same three fixed
starts, at most 200 iterations per start and direction. Keep a point only if
it is inside the domain and its likelihood slack exceeds 1e-6. No truth-guided
repair, retries with new seeds, or selective omission of failed fits.

An attained span exceeding .02 proves that this confidence set cannot support
the old .02 absence-width target for that dataset. A smaller attained span does
not prove a narrow confidence set: numerical searches supply inner witnesses,
not certified global extrema. Excluding old witnesses does not prove absence of
other ambiguous mechanisms. Stop before claiming precision without an outer
enclosure. Record all optimizer outcomes and independently replay saved points
with the repository's high-precision flow arithmetic.

Decision: continue continuous-domain prediction only if the old obstruction is
weakened; prioritize a certified outer envelope next. If >=4/6 new attained
spans still exceed .02, additional information is needed for this procedure.
Neither outcome raises model novelty by itself. Contemporary cofactor/kinetic
models already link microscopic rates with folding; see
[Shah and Tortora et al.](https://doi.org/10.1016/j.molcel.2026.03.037).
