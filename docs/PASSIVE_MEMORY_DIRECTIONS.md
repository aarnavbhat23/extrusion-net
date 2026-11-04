# Extrusion: distinguish stored deformation from ongoing driving

2026-10-30. Spatial minimum-measurement work remains the main priority. This
bounded Extrusion investigation rejects an attractive but insufficient signature
and selects a question worth developing. No new architecture or biological fit
is claimed, and untouched depletion/transfer datasets remain untouched.

## The smallest passive counterexample

After preparation, let two internal deformation modes follow independent passive
Ornstein-Uhlenbeck dynamics

    dx = -A*x*dt + sqrt(2*A)*dW,  A=diag(lambda_1,lambda_2)>0,
    y = x_1+x_2.

Compare initial Gaussian states N((1,-1),I) and N((0,0),I). Both have exactly
the same current visible law y(0)~N(0,2). Their covariance stays I, but their
future visible means differ by

    Delta(t) = exp(-lambda_1*t)-exp(-lambda_2*t).

There is no active driving after preparation. Stored initial deformation alone
produces a future difference. This demonstrates equality of the current visible
law, not equality of the entire pre-preparation history. y is a signed projected
coordinate, not an observed 3D locus-distance distribution or fitted chromatin
model. An isotropic cell population can have zero signed mean despite other
distributional differences; an actual target must handle that geometry.

The algebraic screen evaluates nine rate pairs. All have identical initial
visible laws and nonzero passive future mean differences; peak differences range
from .14815 to .94500 in arbitrary units. A matrix-exponential calculation agrees
with the closed form. These deliberately constructed cases establish the
counterexample, not the frequency or magnitude of memory in cells.

Even a clever contrast is easy to fool. With equally spaced observations and
a=exp(-lambda_1*dt), b=exp(-lambda_2*dt), passive two-mode means obey

    y_2-(a+b)*y_1+a*b*y_0 = 0.

The maximum residual in the nine cases is 1.12e-16. But an unknown constant
measurement offset produces residual offset*(1-a)*(1-b); a third passive mode
also gives a nonzero residual. Both artifacts defeat the two-mode signature in
all nine constructed cases. Rejecting a chosen two-mode model cannot establish
active extrusion. Adding an offset-canceling contrast would also eliminate a
constant equilibrium shift, so extra filtering is not free information.

## A useful conditional boundary

For symmetric passive A with smallest eigenvalue >=lambda_min>0, observation
vector c, common dynamics/noise and each preparation's initial mean norm <=M,
the difference in future mean obeys

    |Delta_y(t)| <= 2*||c||*M*exp(-lambda_min*t).

The proof is Cauchy-Schwarz and ||exp(-At)||_2<=exp(-lambda_min*t). This is a
standard semigroup bound. The assumptions matter: nonnormal dynamics need a
different norm/transient bound; viscoelastic modes, changing mobility and other
active processes are not covered. Unknown M or an arbitrarily slow passive mode
makes a fixed waiting-time argument uninformative. M cannot be bounded from one
current y alone: arbitrarily large opposite internal deformations cancel in y.

One extra calibrated observation z=x_1-x_2 would identify both modes in this toy
model. However, if the rates are already known, y at two different times also
identifies the initial means: the observation matrix [[1,1],[a,b]] has determinant
b-a. The screen verifies this simpler temporal control in every case. It becomes
ill-conditioned when rates are close, and neither claim survives arbitrary
additional hidden modes. This is ordinary observability, not a learning novelty.

## Ranked directions to develop

| Priority | Concrete question | How to try to defeat it | Advance criterion |
|---|---|---|---|
| 1: best scientific question | What additional locus geometry or reference measurement makes a post-depletion response inconsistent with independently bounded passive relaxation? | Add slow Rouse modes, uncertain localization, residual motors, another active force, mobility changes and condition-dependent tracking. Permit every initial state inside the calibrated preparation set. | A useful bound from separate controls; calibrated false attribution; power for a prespecified alternative beyond BILD/Rouse and simple temporal controls; validation by held-out movie/condition. |
| 2: useful observation control | Can crossing elapsed time with illumination dose separate light-driven loss from elapsed-time loss? | Let light alter biology, or add optical dark-state transitions that depend on elapsed time. | An externally measured optical reference and a successful separability check. Treat interval imaging as a control, not a new mechanism. |
| 3: defer architecture claims | Can a compact state represent only distinctions that alter intervention predictions, while retaining uncertainty about unidentifiable motor labels? | Compare predictive-state and state-space models at the same observation/compute budget; require transfer to a truly withheld intervention. | A distinct update rule or guarantee plus capability beyond established predictive-state methods. Renaming a hidden state fails this gate. |

For direction 2, the separable model I(T,D)=S(T)*B(D) implies a 2x2 cross-product
ratio of one. In an exponential specialization, illumination dose can identify
the dose-loss coefficient, while biological and optical elapsed-time losses
remain confounded. Thus even a successful crossed design has a limited target.

Direction 1 has the strongest potential scientific payoff. The possible
contribution would be a useful *minimum-observation requirement for excluding a
specified passive explanation*, with a validated assay. It would not establish
that every rejected passive model implies cohesin activity; attribution to
extrusion additionally needs a selective perturbation and residual-motor bounds.
Nature supplies deformable material with memory, not a license to call all
history dependence a new neural architecture.

## The next bounded experiment

First specify a passive uncertainty class and the exact observable/perturbation.
Use reference conditions only to bound mobility, localization, unresolved modes
and preparation energy/geometry. Check whether one observed locus pair can bound
M at all. If it cannot, identify which second geometry measurement would do so
and whether any public experiment actually supplies it. Two-time extrapolation
with known rates is a mandatory simple control. Do not fit the preparation bound
against the response it is intended to test.

Before evaluating a new condition, freeze movie/cell splits, simultaneous-error
control, a useful effect threshold and a power gate. Include passive models with
extra modes and changed optical offsets as nulls, and an independently specified
ongoing-driving alternative. Compare with the authors' BILD/Rouse inference,
ordinary state-space prediction and a matched temporal control. If the bound is
vacuous, stop this assay as an activity-attribution route and retain the negative
result. If a usable contract is obtained, reserve the independent study for
transfer after acquisition/genomic differences are understood. The present toy
screen is not a preregistration of that future experiment.

## Prior-work limits

[Gabriele et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9069445/) already use
Rouse dynamics in BILD to interpret locus tracking. Their public
[release](https://doi.org/10.5281/zenodo.5770531) has multiple conditions; the
[existing access audit](TRAJECTORY_DATA_FEASIBILITY.md) distinguishes downloaded
baseline data from untouched conditions. Generic latent-loop inference is not new.

[Chan and Rubinstein's theory](https://pmc.ncbi.nlm.nih.gov/articles/PMC10849557/)
already couples extrusion to polymer conformations and relaxation. That linked
version is a preprint; a final publication is available at
[PNAS, DOI 10.1073/pnas.2401494121](https://doi.org/10.1073/pnas.2401494121).
Memory in deformable chromatin is therefore not itself a priority claim.

[Reisser et al.](https://doi.org/10.1038/s41598-020-58634-y) already combine
different dark intervals with controlled illumination to address photobleaching
and dissociation. Crossed time/dose ideas must respect that established protocol.
[Predictive representations of state](https://proceedings.neurips.cc/paper_files/paper/2001/hash/1e4d36177d71bbb3558e43af9577d70e-Abstract.html)
already represent dynamics by action-conditional future predictions. The
counterexample, observability matrix, exponential recurrence and norm bound here
are elementary controls; no original mathematical priority is asserted.

```sh
python study_passive_memory.py
python verify_passive_memory.py
python -m unittest discover -s tests -v
```

All 121 local tests pass. Independent saved-result verification checks the nine
passive examples and 18 offset/extra-mode artifacts without importing the study.

Keep the current assessment near novelty 6.5/10 and research quality 7.5/10,
with roughly one-point judgment uncertainty. This round improves question
selection. It has not supplied biological validation, a new learned mechanism
or grounds to promote the project to 8-9 or submit a manuscript.
