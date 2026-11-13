# Serialization repair

The first frozen run completed all six optimizations, then stopped while writing
JSON because a NumPy boolean was not serializable. Its printed attained spans,
in saved row order, were .04319783462493576, .05636569753071241,
.03234166405258443, .04547848782479358, .021207226960895054,
.04301891923931034. No old witness was retained in any case.

The repair casts the reported boolean and slack to native Python types. Rerun
the identical split seed, inputs, starts, objectives and stopping rules. This
does not constitute a second experiment or fresh confirmation.
