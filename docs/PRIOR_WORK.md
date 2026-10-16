# Prior work and novelty audit

## Biology the mechanism is taken from

- **Loop extrusion.** Sanborn et al., "Chromatin extrusion explains key features of loop and domain formation in wild-type and engineered genomes", PNAS 2015. Fudenberg et al., "Formation of chromosomal domains by loop extrusion", Cell Reports 2016.
- **Direct observation of cohesin extruding DNA.** Davidson et al., "DNA loop extrusion by human cohesin", Science 2019. Kim et al., "Human cohesin compacts DNA by loop extrusion", Science 2019.
- **Convergent CTCF orientation.** Rao et al., "A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping", Cell 2014. de Wit et al., "CTCF binding polarity determines chromatin looping", Molecular Cell 2015. Guo et al., "CRISPR inversion of CTCF sites alters genome topology and enhancer/promoter function", Cell 2015.

The model's stall rule (forward sites stop left-moving arms, reverse sites stop right-moving arms) and the task's insulation rule both follow this work. The model simplifies the biology heavily: one motor at a time, independent arms, no collisions, no unloading beyond stalling, and no 3D polymer physics.

## Closest machine-learning work

| Work | Relation | Difference |
|---|---|---|
| dLEM, [bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.10.17.682904.full.pdf) | Makes loop extrusion differentiable, with position-specific left/right velocity profiles fitted to Hi-C. | dLEM predicts genome folding. Here extrusion is the information-routing step of a general sequence model, trained only on downstream labels and never shown a contact map. |
| Stick-breaking attention, Tan et al., ICLR 2025 ([arXiv 2410.17980](https://arxiv.org/abs/2410.17980)) | Attention weights from cumulative products of stopping probabilities. | Its break probabilities depend on query-key pairs and are one-sided from the query. Extrusion uses per-position stalls, two independent arms from a separate load site, and a sum over load sites. The `one_arm` ablation is the nearest form to stick-breaking in this repo. |
| Forgetting Transformer, Lin et al., ICLR 2025 | Cumulative forget gates added to attention logits. | Same family of cumulative-product decays; no load sites or oriented barriers. |
| ALiBi, Press et al., ICLR 2022 ([arXiv 2108.12409](https://arxiv.org/abs/2108.12409)); RoPE, Su et al. ([arXiv 2104.09864](https://arxiv.org/abs/2104.09864)) | Position schemes for attention, ALiBi built for length extrapolation. | Used here as controls. |
| Sequence-to-Hi-C predictors, e.g. Akita (Fudenberg, Kelley and Pollard, Nature Methods 2020) | Deep models that predict contact maps from DNA sequence. | They predict chromatin structure with generic layers. This repo uses an extrusion process as a layer. |

## Possible overlap not yet checked

- **Span-based parsing and semi-Markov models.** `D(i, j)` sums, over split points `k`, products of left and right quantities. That is the shape of inside-algorithm and segment scores in parsing. A dedicated search for "differentiable span insulation" or "barrier-gated span attention" has not been done.
- **Gated linear attention and selective state-space models.** Cumulative gates are standard there; a two-sided, load-site version may exist under another name.
- **Conference proceedings.** Only web searches were run (October 2026). NeurIPS, ICML and ICLR proceedings and the bioRxiv genomics-ML literature were not searched systematically.

## What would count as a first

"First" is not claimed. The defensible statement is narrower: this is an implementation and controlled test of two-arm, oriented-barrier loop extrusion as a sequence model's information-routing operator, with a released ablation against its nearest one-arm form.
