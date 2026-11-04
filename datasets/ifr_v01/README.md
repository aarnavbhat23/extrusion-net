# Reproduce the measured iFRAP curation

This small source-only package reproduces the Extrusion portion of the local
measurement dataset release. It downloads the original 4.52 MB public CSV,
enforces its publisher SHA-256, preserves attribution, flags unresolved source
issues and builds the same whole-day partition. Large/generated data stay ignored.
The processing references are pinned to their source commit and retain the
BSD 3-Clause license in `source/ifr_LICENSE`.

From the repository root, after installing the repository's requirements:

```sh
python datasets/ifr_v01/build_extrusion.py
python study_measured_ifrap.py --dataset-root datasets/ifr_v01
python verify_measured_ifrap.py --dataset-root datasets/ifr_v01
```

The [dataset card](extrusion/DATASET_CARD.md) documents observable limitations,
the eight quarantined movies and absent biological ground truth. The frozen
[study protocol](../../docs/MEASURED_IFRAP_PROTOCOL.md) controls evaluation.
The source is [Louafi and Giorgetti, version 1](https://doi.org/10.17632/spp4x8jpdy.1),
CC BY 4.0; these are derived public measurements, not newly collected experiments.
