# Deterministic four-response screen

Exploratory theory-led protocol written before running this experiment. Retain
the earlier random-probe screens. Use the same two complementary alternating
loads q0=.5, contrast .8, and two value vectors: all ones and the repeated code
`++++----`. This separates every distance-two edge derivative at zero gates.
Derive a length-independent operator-norm tail bound for gates in [0,rho]. A
positive lower Lipschitz constant would imply global, rather than only local,
identification on that restricted domain.

Check the Jacobian and global secant bounds at 8/16/32/64 loads, 10 mechanisms
per setting, rho=0/.02/.05/.08/.2/.8. Compare random sign value probes at the
same loading and response budgets. Include the regimes where the proved bound
is negative; they do not inherit the theorem. Do not infer necessity from a
sufficient condition.

Fit free gates from four noisy vector responses at 8/16/32 loads, five mechanisms
per setting, rho=.05, Gaussian noise SD .001/.01, three fixed starts and matched
random/deterministic probes. Record actual recovery, residual, held-out loading
responses, convergence, and the computable global recovery-error bound from a
95% Gaussian noise-norm bound. A candidate's residual gives a bound regardless
of whether a nonlinear optimizer is globally optimal. Input loads are known,
transport is unchanged across probes, and output noises are independent.

Seed 20261108 is a reproducibility identifier, not the study date. This is a
specialized recovery theorem and matched-generator screen. General machine
learning advantages and priority require further evidence.
