# Equation and protocol originality audit — 2026-10-30

This is a targeted comparison, not an exhaustive priority review. The proposed
counterexample was specified before this audit; finding familiar mathematical
structure is a reason to narrow the claim, not rename the same idea.

| Closest work | Equation/protocol overlap | What this round actually adds |
|---|---|---|
| [Darwin and Kiefer, HMM equivalence, 2020](https://doi.org/10.4230/LIPIcs.FSTTCS.2020.43) | Equality of observation laws is the core object; their general equivalence machinery is far beyond our identical-factor example. | A specific controlled biological/optical interpretation counterexample and its reporter assumptions. No new equivalence algorithm or theorem. |
| [Bryan et al., binding kinetics and photophysics, 2025](https://doi.org/10.1021/acs.jpcb.5c01180) | Joint consideration of binding/unbinding and photobleaching; fluorescence intensity is not a direct binding-state label. Their method infers rates with an HMM-based pipeline and includes DNA-origami experiments. | An adversarial unrestricted optical-response boundary, not evidence that their constrained generative model is unidentifiable. We have not implemented or beaten their estimator. |
| [Minimally sufficient experimental design using identifiability analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC10771435/) | Structural analysis followed by profile-likelihood practical identification and measurement selection; weighted residual objectives are conventional. | A restricted calibrated contrast and a direct-reporter control. It does not replace their general workflow. |
| [Existing repository pulse-order study](PULSE_ORDER_RESULTS.md) | Balanced schedules, scalar-hazard cancellation, hidden-state alternatives, corruption bounds and failed power tests already exist here. | Attribution to biology versus optical hidden dynamics, with an explicit extra-reporter failure. Pulse order itself is not a new direction. |

The confidence radius is Hoeffding's inequality; interval inversion is standard
set identification. The all-schedule equivalence follows immediately from equal
transition/emission factors, including under an adaptive policy. Neither should
be presented as a novel mathematical tool.

## Decision

The generic 'order reveals biological memory' claim is not defensible. A useful
research opening would be an experimentally justified restriction that excludes
optical equivalence and yields calibrated prediction under unknown intervention
strength. To claim originality, specify how that restriction and resulting
observable differ from binding/photobleaching inference and conventional optimal
design, implement the closest compatible estimator, and obtain independent
state validation. No priority or 8–9 novelty rating is established this round.
