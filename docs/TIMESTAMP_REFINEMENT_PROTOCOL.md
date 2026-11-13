# Stronger timestamp control — frozen 2026-11-01

The first loading-target screen used separate extrema for L=log(ell) and
S=A+W. That outer interval is valid but discards their dependence. Its failed
precision gate cannot establish that timestamp data are less useful. Freeze this
development refinement before evaluating its outcomes. Keep the original result.

Use all 300 previously saved matched N=24000 timestamp experiments, with their
same exact binomial rectangles and f in [1.8,2.2]. This reuses observed development
data and is not a fresh confirmation experiment or a change to the original gate.

Find min/max feasible L by the existing linear constraints. Divide that full
interval into exactly 16 equal closed bins (or one if endpoints agree). In every
bin solve min/max S. Its target presence lies in

    [exp(L_bin_low)/(1+S_max), exp(L_bin_high)/(1+S_min)].

Union the bin intervals, then transform to absence. Skip solver-infeasible bins;
if no bin remains use conservative [0,1] and report status. Save all LP primal/
dual witnesses and bin limits. This is a certified outer bound, not a sharp
timestamp optimum. For bin width w, each one-sided probability slack from this
factorization is bounded by 1-exp(-w), since presence<=1. Report that numerical
relaxation tolerance, coverage and width for each mechanism, and compare to the
snapshot methods on the same resource ledger. No model or acquisition-time
choices are retuned. A strong control winning is retained as the result.
