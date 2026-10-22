# Input-integrity amendment before an evaluable primary result

The first frozen run stopped on duplicate (movie, track, frame) keys; it produced
no p-value or result artifact. A metadata-only key audit found 9 repeated keys
in 2 unfiltered tracks, and 1 repeated key in 1 filtered track. All frame indices
are integers. Choosing one distance among duplicate localizations would introduce
an undocumented tracking decision.

Exclude each entire affected track in both pipelines, recording its identifier
and row count. This rule depends only on observation keys and is equivariant
under time reversal. Keep the frozen lag scores, split, thresholds, cluster unit,
randomization and significance gate. This is an explicitly amended screen, not an
unchanged preregistration. The excluded tracks and unidentified acquisition dates
remain a limitation, alongside the unverified independence of acquisition dates.
