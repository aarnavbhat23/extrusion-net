"""Exact-arithmetic audit of the counted-reference counterexample and the
loading-drift sharpness proof (docs/LOADING_TARGET_THEORY.md,
docs/COUNTED_LOADING_THEORY.md at d9cd4d3).

Self-contained: imports nothing from the repository, uses no optimizer and no
floating point in any check that is reported as exact. Run with any Python 3.
"""
from fractions import Fraction as F
from itertools import combinations
import json
import random

GAMMA, TARGET = F(3, 10), F(3)


# ---------------------------------------------------------------- the model
def presence(ell, mu, delta, f, gamma=GAMMA, kappa=F(1)):
    """Absolute presence at an Exp(gamma) clock, h=1, loading kappa*ell."""
    return kappa*ell*gamma/(gamma+delta+f*mu)


def target_absence(ell, mu, delta, g=TARGET, gamma=GAMMA):
    return 1-presence(ell, mu, delta, g, gamma)


# ------------------------------------------------- A. the 9/286 counterexample
def audit_counterexample():
    ell, f, q0 = F(3, 4), F(2), F(9, 10)
    worlds = dict(A=dict(mu=F(1, 5), delta=F(1, 5), kappa=F(4, 5), q1=F(9, 10)),
                  B=dict(mu=F(3, 10), delta=F(1, 10), kappa=F(8, 9), q1=F(81, 100)))
    out = {}
    for name, w in worlds.items():
        assert w['mu'] >= 0 and w['delta'] >= 0 and 0 < w['kappa'] <= 1 and 0 < w['q1'] <= 1
        obs = [presence(ell, w['mu'], w['delta'], 1),
               presence(ell, w['mu'], w['delta'], f, kappa=w['kappa']),
               q0*ell, w['q1']*w['kappa']*ell]
        out[name] = dict(observed=obs, target=target_absence(ell, w['mu'], w['delta']),
                         sensitivity_ratio=w['q1']/q0)
    assert out['A']['observed'] == out['B']['observed'] == [F(9, 28), F(1, 5), F(27, 40), F(27, 50)]
    gap = out['B']['target']-out['A']['target']
    assert gap == F(9, 286)

    # The two worlds are exactly the endpoints of the exact-data identified set
    # when f=2 is known and q1/q0 in [9/10,1]: kappa=(B1/B0)/r, J=kappa/R1,
    # W=(g-1)(J-U)/(f-1), target presence 1/(U+W).
    R0, R1, B0, B1 = out['A']['observed']
    U = 1/R0
    ends = []
    for r in (F(1), F(9, 10)):
        J = (B1/B0)/r/R1
        W = (TARGET-1)*(J-U)/(f-1)
        assert W <= (TARGET-1)*(U-1)  # hazard constraint, ell_max=1
        ends.append(1-1/(U+W))
    assert ends == [out['A']['target'], out['B']['target']]

    # Dwell-time readout distinguishes them: calibration-arm loss hazards differ.
    hazards = {n: worlds[n]['delta']+f*worlds[n]['mu'] for n in worlds}
    assert hazards['A'] != hazards['B']
    return dict(observed=[str(x) for x in out['A']['observed']],
                targets={n: str(out[n]['target']) for n in out}, gap=str(gap), gap_pp=float(gap)*100,
                worlds_are_endpoints_of_exact_set_f_known_ratio_09_to_1=True,
                calibration_arm_hazards={n: str(h) for n, h in hazards.items()})


# ------------------------------- B. necessity and sufficiency of the polytope
def polytope_ok(U, V, J, W, kL, kU, fL, fU, ell_max, g=TARGET):
    tL, tU = (g-1)/(fU-1), (g-1)/(fL-1)
    return (kL*V <= J <= kU*V and J >= U and tL*(J-U) <= W <= tU*(J-U)
            and U >= 1/ell_max and W <= (g-1)*(U-1/ell_max) and W >= 0)


def forward(ell, mu, delta, kappa, f, gamma=GAMMA, g=TARGET):
    U = 1/presence(ell, mu, delta, 1, gamma)
    V = 1/presence(ell, mu, delta, f, gamma, kappa)
    return U, V, kappa*V, (g-1)*mu/(ell*gamma)


def reconstruct(U, V, J, W, ell_max, gamma=GAMMA, g=TARGET):
    mu = ell_max*gamma*W/(g-1)
    delta = gamma*(ell_max*U-1)-mu
    kappa = J/V
    f = None if W == 0 else 1+(g-1)*(J-U)/W
    return mu, delta, kappa, f


def rand_frac(rng, lo, hi, den=997):
    return lo+(hi-lo)*F(rng.randrange(den+1), den)


def audit_sharpness(trials=4000, seed=20261103):
    rng = random.Random(seed)
    kL, kU, fL, fU = F(4, 5), F(1), F(9, 5), F(11, 5)
    n_forward = n_reverse = n_boundary = 0
    for _ in range(trials):
        ell_max = rand_frac(rng, F(1, 2), F(1))
        ell = rand_frac(rng, F(1, 20), ell_max)
        mu, delta = rand_frac(rng, F(0), F(2)), rand_frac(rng, F(0), F(2))
        kappa, f = rand_frac(rng, kL, kU), rand_frac(rng, fL, fU)
        U, V, J, W = forward(ell, mu, delta, kappa, f)
        assert polytope_ok(U, V, J, W, kL, kU, fL, fU, ell_max)          # necessity
        assert 1/(U+W) == presence(ell, mu, delta, TARGET)              # target identity
        n_forward += 1
    # Sufficiency: sample feasible points directly (including boundary faces).
    while n_reverse < trials:
        ell_max = rand_frac(rng, F(1, 2), F(1))
        U = 1/ell_max+rand_frac(rng, F(0), F(4))
        V = rand_frac(rng, U, U/kL+F(4))
        lo, hi = max(U, kL*V), kU*V
        if lo > hi:
            continue
        jface, wface = rng.randrange(3), rng.randrange(3)
        J = lo if jface == 0 else hi if jface == 1 else rand_frac(rng, lo, hi)
        tL, tU = (TARGET-1)/(fU-1), (TARGET-1)/(fL-1)
        wlo, whi = tL*(J-U), min(tU*(J-U), (TARGET-1)*(U-1/ell_max))
        if wlo > whi:
            continue
        W = wlo if wface == 0 else whi if wface == 1 else rand_frac(rng, wlo, whi)
        assert polytope_ok(U, V, J, W, kL, kU, fL, fU, ell_max)
        mu, delta, kappa, f = reconstruct(U, V, J, W, ell_max)
        assert mu >= 0 and delta >= 0 and kL <= kappa <= kU and kappa*ell_max <= 1
        if W == 0:
            assert J == U and mu == 0
            f = fL  # irrelevant when mu=0
            n_boundary += 1
        assert fL <= f <= fU
        assert 1/presence(ell_max, mu, delta, 1) == U
        assert 1/presence(ell_max, mu, delta, f, kappa=kappa) == V
        assert 1/presence(ell_max, mu, delta, TARGET) == U+W
        n_reverse += 1
    return dict(forward_mechanisms_checked=n_forward, feasible_points_reconstructed=n_reverse,
                W_zero_boundary_points=n_boundary, exact_rational=True)


# --------------------- C. exact LP extremes by vertex enumeration (no solver)
def solve(A, b):
    n = len(A)
    M = [list(r)+[v] for r, v in zip(A, b)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                k = M[r][c]/M[c][c]
                M[r] = [x-k*y for x, y in zip(M[r], M[c])]
    return [M[i][n]/M[i][i] for i in range(n)]


def drift_lp_exact(R0box, R1box, kL, kU, fL, fU, ell_max, g=TARGET):
    """All constraints as a.x <= b on x=(U,V,J,W); enumerate vertices."""
    tL, tU = (g-1)/(fU-1), (g-1)/(fL-1)
    rows = [([0, kL, -1, 0], 0), ([0, -kU, 1, 0], 0), ([1, 0, -1, 0], 0),
            ([-tL, 0, tL, -1], 0), ([tU, 0, -tU, 1], 0), ([-(g-1), 0, 0, 1], -(g-1)/ell_max),
            ([-1, 0, 0, 0], -max(1/R0box[1], 1/ell_max)), ([1, 0, 0, 0], 1/R0box[0]),
            ([0, -1, 0, 0], -1/R1box[1]), ([0, 1, 0, 0], 1/R1box[0]),
            ([0, 0, -1, 0], 0), ([0, 0, 0, -1], 0)]
    rows = [([F(x) for x in a], F(b)) for a, b in rows]
    values = []
    for combo in combinations(rows, 4):
        x = solve([a for a, _ in combo], [b for _, b in combo])
        if x is None:
            continue
        if all(sum(ai*xi for ai, xi in zip(a, x)) <= b for a, b in rows):
            values.append(x[0]+x[3])
    if not values:
        return None
    return 1-1/min(values), 1-1/max(values)


def closed_form_exact(R0, R1, kL, kU, fL, fU, ell_max, g=TARGET):
    """Exact-data (degenerate rectangle) identified interval, closed form."""
    U, V = 1/R0, 1/R1
    tL, tU = (g-1)/(fU-1), (g-1)/(fL-1)
    Jmin, Jmax = max(U, kL*V), kU*V
    if Jmin > Jmax or U < 1/ell_max:
        return None
    cap = (g-1)*(U-1/ell_max)
    Wmin = tL*(Jmin-U)
    if Wmin > cap:
        return None
    Wmax = min(tU*(Jmax-U), cap)
    return 1-1/(U+Wmin), 1-1/(U+Wmax)


def audit_lp(seed=7):
    rng = random.Random(seed)
    checks = 0
    for _ in range(40):
        ell = rand_frac(rng, F(1, 2), F(1), 50)
        mu, delta = rand_frac(rng, F(1, 20), F(1), 50), rand_frac(rng, F(1, 20), F(1), 50)
        kappa, f = rand_frac(rng, F(4, 5), F(1), 50), rand_frac(rng, F(9, 5), F(11, 5), 50)
        R0, R1 = presence(ell, mu, delta, 1), presence(ell, mu, delta, f, kappa=kappa)
        exact = closed_form_exact(R0, R1, F(4, 5), F(1), F(9, 5), F(11, 5), F(1))
        lp = drift_lp_exact((R0, R0), (R1, R1), F(4, 5), F(1), F(9, 5), F(11, 5), F(1))
        assert exact == lp, (exact, lp)
        truth = target_absence(ell, mu, delta)
        assert exact[0] <= truth <= exact[1]
        checks += 1
    return dict(exact_data_cases=checks, closed_form_equals_vertex_enumeration=True)


# ---------------------------------- D. boundaries outside the declared scope
def audit_false_positive():
    """Equal sensitivity, but a reference false-positive rate phi among
    unloaded episodes. Matches world A's four probabilities exactly."""
    ell, phi = F(3, 4), F(1, 20)
    q = (F(27, 40)-phi*(1-ell))/ell
    # B1 = q*kappa*ell + phi*(1-kappa*ell) = 27/50  ->  solve for kappa
    kappa = (F(27, 50)-phi)/(q*ell-phi*ell)
    mu = F(9, 8)*kappa-F(7, 10)
    delta = F(2, 5)-mu
    assert 0 <= mu and 0 <= delta and 0 < kappa <= 1 and 0 < q <= 1
    obs = [presence(ell, mu, delta, 1), presence(ell, mu, delta, 2, kappa=kappa),
           q*ell+phi*(1-ell), q*kappa*ell+phi*(1-kappa*ell)]
    assert obs == [F(9, 28), F(1, 5), F(27, 40), F(27, 50)]
    shift = target_absence(F(3, 4), F(1, 5), F(1, 5))-target_absence(ell, mu, delta)
    return dict(phi=str(phi), q=str(q), kappa=str(kappa), target=str(target_absence(ell, mu, delta)),
                shift_vs_world_A_pp=float(shift)*100,
                note='naive ratio B1/B0 reads 4/5; true kappa is 98/125')


def audit_kappa_above_one():
    """The proof only uses kappa*ell_max<=1, not kappa<=1."""
    ell_max, kU = F(4, 5), F(5, 4)
    point = (F(3), F(3), F(15, 4), F(3, 2))
    assert polytope_ok(*point, F(4, 5), kU, F(9, 5), F(11, 5), ell_max)
    mu, delta, kappa, f = reconstruct(*point, ell_max)
    assert kappa == kU and kappa*ell_max == 1 and mu >= 0 and delta >= 0 and f == 2
    assert 1/presence(ell_max, mu, delta, f, kappa=kappa) == point[1]
    # With ell_max=1 the same point would need treated loading 5/4: not physical.
    return dict(condition_used_by_proof='kU*ell_max <= 1',
                example=dict(ell_max='4/5', kU='5/4', treated_loading=str(kappa*ell_max)),
                when_violated='LP remains a valid outer image (every constraint is still necessary) '
                              'but the sharp set gains the bilinear constraint U*V >= J')


def audit_state_resolved():
    """Does adding a retained motor-state readout break the counterexample?
    Two loaded states, start in state 1, shared generator across arms.
    Normalized share of state 2 at total exit rate c: p/(c+p+r)."""
    c0, cA, cB = F(7, 10), F(9, 10), F(1)   # gamma+a in arm 0 (both), arm 1 (A, B)
    share = lambda p, r, c: p/(c+p+r)
    # (i) unknown reversible dynamics: an exact matching world exists.
    pA, rA = F(1, 2), F(1, 5)
    pB, rB = pA*(cB-c0)/(cA-c0), (cB-c0)/(cA-c0)*(rA+c0)-c0
    assert rB >= 0
    assert share(pA, rA, c0) == share(pB, rB, c0) and share(pA, rA, cA) == share(pB, rB, cB)
    # (ii) irreversible two-state (r=0): the share pins p/c in each arm, so
    # c1/c0 -- and with R1/R0 the loading ratio -- is identified.
    sA0, sA1 = share(pA, 0, c0), share(pA, 0, cA)
    ratio = (sA0/(1-sA0))/(sA1/(1-sA1))   # = c1/c0
    assert ratio == cA/c0
    # kappa = (R1/R0)*(c1/c0) recovers world A's 4/5 without any reference.
    kappa = (F(1, 5)/F(9, 28))*ratio
    assert kappa == F(4, 5)
    return dict(reversible_unknown_generator='counterexample survives (exact matching world)',
                world_B_generator=dict(p=str(pB), r=str(rB)),
                irreversible_or_calibrated_generator='loading ratio identified internally; kappa=4/5 recovered',
                caveat='requires the perturbation to leave state dynamics unchanged')


def structural_floor():
    """Exact-data widths for the frozen counted-study mechanisms: which
    contract dominates the width before any sampling error?"""
    mechs = [(F(1, 5), F(1, 5), F(3, 4)), (F(2, 5), F(1, 10), F(17, 20)), (F(13, 20), F(3, 25), F(19, 20))]
    rows = []
    for (mu, delta, ell) in mechs:
        for kappa in (F(4, 5), F(9, 10), F(1)):
            R0, R1 = presence(ell, mu, delta, 1), presence(ell, mu, delta, 2, kappa=kappa)
            def width(fL, fU, kL, kU):
                lo, hi = closed_form_exact(R0, R1, kL, kU, fL, fU, F(1))
                return float(hi-lo)*100
            kq = (kappa/F(1), min(F(1), kappa/F(9, 10)))   # ratio contract [0.9,1]
            rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], kappa=str(kappa),
                             target_absence=float(target_absence(ell, mu, delta)),
                             f_uncertain_kappa_exact=width(F(9, 5), F(11, 5), kappa, kappa),
                             f_exact_sensitivity_ratio_uncertain=width(F(2), F(2), *kq),
                             both=width(F(9, 5), F(11, 5), *kq),
                             f_within_2pct_sensitivity_exact=width(F(49, 25), F(51, 25), kappa, kappa)))
    return rows


def naive_extrapolation_gaps():
    """Distance from the true target to two naive extrapolations (kappa=1):
    'no extra effect beyond f' and 'linear in dose on the absence scale'."""
    mechs = [(F(1, 5), F(1, 5), F(3, 4)), (F(2, 5), F(1, 10), F(17, 20)), (F(13, 20), F(3, 25), F(19, 20))]
    rows = []
    for mu, delta, ell in mechs:
        a0, a1, ag = (1-presence(ell, mu, delta, f) for f in (1, 2, 3))
        linear = a1+(a1-a0)  # g=3 is one more dose step beyond f=2
        rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], baseline=float(a0),
                         calibration_arm=float(a1), target=float(ag),
                         gap_to_no_extra_effect_pp=float(ag-a1)*100,
                         gap_to_linear_extrapolation_pp=float(linear-ag)*100))
    return rows


MECHS = [(F(1, 5), F(1, 5), F(3, 4)), (F(2, 5), F(1, 10), F(17, 20)), (F(13, 20), F(3, 25), F(19, 20))]


def exact_width_general(R0, R1, fL, fU, gL, gU, kL, kU, ell_max):
    """Exact-data identified interval when the TARGET multiplier g is also
    only known to [gL,gU]. M=mu/(ell*gamma): J-U=(f-1)M, target 1/R=U+(g-1)M."""
    U, V = 1/R0, 1/R1
    Jmin, Jmax = max(U, kL*V), kU*V
    cap = U-1/ell_max
    Mmin = (Jmin-U)/(fU-1)
    Mmax = min((Jmax-U)/(fL-1), cap)
    assert Jmin <= Jmax and Mmin <= cap
    lo, hi = 1-1/(U+(gL-1)*Mmin), 1-1/(U+(gU-1)*Mmax)
    return float(hi-lo)*100


def tolerance_scan():
    """Exact-data width (pp) vs relative calibration tolerance e on the factor
    (f in 2(1+-e)), with the target multiplier g exact or equally uncertain,
    and vs the sensitivity-ratio lower bound r (q1/q0 in [r,1], true ratio 1)."""
    rows = []
    for mu, delta, ell in MECHS:
        kappa = F(4, 5)
        R0, R1 = presence(ell, mu, delta, 1), presence(ell, mu, delta, 2, kappa=kappa)
        for pct in (0, 1, 2, 4, 6, 10):
            e = F(pct, 100)
            rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], kind='factor', tolerance_pct=pct,
                             g_exact=exact_width_general(R0, R1, 2*(1-e), 2*(1+e), F(3), F(3), kappa, kappa, F(1)),
                             g_same_relative_error=exact_width_general(R0, R1, 2*(1-e), 2*(1+e), 3*(1-e), 3*(1+e), kappa, kappa, F(1))))
        for r in (F(1), F(99, 100), F(97, 100), F(95, 100), F(9, 10)):
            rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], kind='sensitivity_ratio_lower', r=float(r),
                             width_f_exact=exact_width_general(R0, R1, F(2), F(2), F(3), F(3), kappa, min(F(1), kappa/r), F(1))))
    return rows


def kappa_cap_dependence():
    """At true kappa=1 the sensitivity contract collapses only because kappa<=1
    is imposed. Compare cap 1 with cap 1/ell_max (ell_max=19/20, still sharp)."""
    rows = []
    ell_max = F(19, 20)
    for mu, delta, ell in MECHS:
        R0, R1 = presence(ell, mu, delta, 1), presence(ell, mu, delta, 2)
        for cap in (F(1), 1/ell_max):
            rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], kappa_cap=float(cap),
                             width_f_exact_ratio_in_09_1=exact_width_general(R0, R1, F(2), F(2), F(3), F(3), F(1), min(cap, F(10, 9)), ell_max)))
    return rows


def combined_tolerance():
    """True kappa=1; f and g share relative error e; arm loading ratio in
    [1-l, 1+l] with ell_max=1/(1+l) so that kU*ell_max<=1 keeps the LP sharp."""
    rows = []
    for mu, delta, ell in MECHS:
        R0, R1 = presence(ell, mu, delta, 1), presence(ell, mu, delta, 2)
        for pf, pl in [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3), (2, 5), (4, 3), (4, 5)]:
            e, l = F(pf, 100), F(pl, 100)
            rows.append(dict(mechanism=[str(mu), str(delta), str(ell)], fg_tolerance_pct=pf, loading_tolerance_pct=pl,
                             width=exact_width_general(R0, R1, 2*(1-e), 2*(1+e), 3*(1-e), 3*(1+e), 1-l, 1+l, 1/(1+l))))
    return rows


if __name__ == '__main__':
    report = dict(counterexample=audit_counterexample(), sharpness=audit_sharpness(),
                  combined_tolerance_pp=combined_tolerance(),
                  tolerance_scan_pp=tolerance_scan(), kappa_cap_dependence_pp=kappa_cap_dependence(),
                  lp_vertex_check=audit_lp(), false_positive_boundary=audit_false_positive(),
                  kappa_above_one=audit_kappa_above_one(), state_resolved=audit_state_resolved(),
                  structural_width_floor_pp=structural_floor(),
                  naive_extrapolation_gaps=naive_extrapolation_gaps())
    print(json.dumps(report, indent=2))
