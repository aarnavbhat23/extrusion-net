# Loading-pattern screen

Exploratory protocol, written before this screen's results. The question is
whether changing the source pattern, while holding transport gates fixed, removes
directional ambiguity more efficiently than local gate interventions. This is an
input-design experiment, not a new neural architecture.

Compare two repeated uniform patterns, two random complementary binary patterns,
and two alternating complementary patterns. All use the same number of contact
observations and the same total loading over the pair. Alternating and random
patterns also have identical squared loading energy. Uniform has lower squared
energy; equal total loading is not equal energy.

1. Prove and check local two-view inversion against the original contact sum.
   Use 200 mechanisms at each of 8, 32 and 128 loads, gates uniform in [.2,.8],
   known loads centered at .5 with contrast .8, Gaussian contact noise SD .01.
   Report unresolved gates, untruncated gate MSE, and failure when gates change
   between views. Negative/inaccurate recovered gates remain in the score.
2. Test a shared four-feature logistic controller, rather than extending the
   free-gate result to it. Ten seeds, 16 loads, one training feature context,
   independent length-64 test context, both nearest-band and full-map fitting.
   Match two views and noise across patterns; use three fixed optimization starts.
   Report gate prediction and local right-cut response errors. Full-map repeated
   uniform input is a necessary strong control. No task-label learning claim.

Seed 20261105 is an arbitrary reproducibility identifier, not a study date.
Continue only if the exact result is correct and the shared/full-map control
leaves a useful effect. Even a positive screen needs externally observable probe
responses; a model cannot supervise itself by inspecting its own hidden gates.
If this requires unavailable supervision, retain the result as operator theory
and pivot the learning hypothesis to experimentally testable routing transfer.

## Follow-up screen after the contact-observation pilot

Before running the following check: replace full contact-map observations with
two or four vector responses `D(q) v`. This asks whether a linear number of
observable responses can distinguish the free gates. Use 10 mechanisms at sizes
8/16/32/64, gates in [.2,.8], repeated uniform versus alternating loads, one or
two independently sampled Rademacher value probes per pattern. Match the value
probes across loading strategies. Report local Jacobian rank, smallest singular
value and mean local inverse-information variance, not just rank. Reject a claim
of stable identification if variance grows badly with length. This is a local
screen; it neither proves global identification nor demonstrates training.

After the local screen, test actual fitting from four vector responses (68 noisy
scalar observations) with the same four-feature shared logistic rule. Ten new
synthetic seeds within the exploratory stream, train length 16 and test length
64; uniform, random-complementary and alternating loading, identical probes and
noise across strategies, three fixed starts. Score the unseen-context gates and
local right-cut responses. This still requires an externally observed system;
it is not supervision from the network's own predictions.

After fitting the vector-response controller, the random loading control beat
the alternating code on average. Add random complementary loading to the local
information screen as well, using a separate RNG stream so existing mechanisms,
probes and shared-controller fits remain unchanged. This tests whether source
diversity, rather than the specific alternating pattern, explains the result.
