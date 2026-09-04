# RATIO-DIA

**R**esponse-**A**ssociated **T**opology **I**ntegration and **O**rganization for
data-independent acquisition mass spectrometry.

RATIO-DIA is a reproducible workflow for converting an overconnected DIA-MS/MS
molecular network into topology-focused spectral modules and prioritizing
features associated with a bioactive fraction. It was developed for complex
natural-product mixtures in which co-fragmentation and shared diagnostic ions
can produce a giant connected component.

## What the workflow does

1. Reads centroided MGF spectra.
2. Calculates modified-cosine similarity or reuses an existing edge table.
3. Applies a reciprocal Top-*K* neighborhood constraint.
4. Detects non-overlapping topology-focused spectral modules.
5. Calculates fraction-level recovery indices from biological responses.
6. Integrates feature abundance, target-fraction enrichment, and agreement with
   the activity gradient to assign Tier 1 and Tier 2 candidates.
7. Exports node, edge, module-summary, and focused-network tables for Cytoscape.
8. Optionally evaluates module structural coherence using all CHON candidate
   structures rather than only the top-ranked annotation.

## Installation

```bash
git clone https://github.com/christalho22/RATIO-DIA.git
cd RATIO-DIA
python -m venv .venv
```

Activate the environment and install:

```bash
pip install -e .
```

For candidate-structure consistency analysis:

```bash
pip install -e ".[structure]"
```

## Quick start

Run the included miniature example:

```bash
ratio-dia \
  --mgf examples/demo.mgf \
  --abundance-table examples/abundance_table.csv \
  --activity-table examples/activity_response.csv \
  --output outputs/demo \
  --min-cosine 0.70 \
  --min-matched-peaks 6 \
  --top-k 5 \
  --top-k-mode mutual \
  --mode weak \
  --hcpin-threshold 1.0 \
  --min-module-size 3
```

Windows PowerShell uses backticks instead of backslashes for multiline commands,
or the command can be entered on one line.

For a large study, calculating every spectral pair is quadratic. Reuse a
validated similarity table whenever possible:

```bash
ratio-dia --mgf study.mgf --edges-in modified_cosine_edges.csv \
  --abundance-table abundance.csv --activity-table activity.csv \
  --output outputs/study
```

## Main outputs

| File | Purpose |
|---|---|
| `01_filtered_similarity_edges.csv` | Similarity-thresholded network before Top-*K* focusing |
| `02_nodes_with_modules_and_D_labels.csv` | Node attributes, module IDs, and activity scores |
| `03_edges_with_modules.csv` | Focused edges with topology and module assignments |
| `04_module_summary.csv` | Module size, density, mean similarity, and fraction abundance |
| `05_D_focused_nodes.csv` | Nodes in modules containing Tier 1/2 candidates |
| `06_D_focused_edges.csv` | Internal edges of candidate-containing modules |
| `07_D_active_marker_mz_by_module.csv` | Prioritized candidate features by module |

For Cytoscape, import `03_edges_with_modules.csv` as an undirected network using
`Source_Scan` and `Target_Scan`, then import `02_nodes_with_modules_and_D_labels.csv`
as node-table columns using `Alignment_ID` as the key. Color nodes by `Module_ID`
and highlight `D_candidate = 1`.

## Structural-coherence analysis

Prepare a long candidate table using
`examples/candidate_structures_template.csv`, then run:

```bash
ratio-dia-structure \
  --nodes outputs/study/02_nodes_with_modules_and_D_labels.csv \
  --candidates candidate_structures.csv \
  --min-module-size 11 \
  --output outputs/study/module_structure_consistency.csv
```

The analysis reports annotation coverage, maximum all-candidate Morgan-fingerprint
Tanimoto similarity, shared Murcko scaffolds, and shared candidate ontologies.
An "any candidate" match is evidence of structural interpretability, not a
confirmed identification.

## Reproducibility and interpretation

- Record the MGF export settings, mass tolerances, similarity threshold, minimum
  matched peaks, Top-*K*, module condition, and minimum module size.
- `Module_ID = 0` denotes nodes not assigned to a module meeting the minimum size.
- A topology-focused spectral module is not automatically a confirmed compound
  family or a causal bioactivity mechanism.
- DIA co-fragmentation, adducts, in-source fragments, non-covalent aggregates,
  and co-eluting isomers can create redundant or ambiguous nodes.
- Fraction-level correlations are intended for prioritization. Isolated compounds
  require orthogonal validation by authentic standards, HRMS/MS, NMR, and
  biological assays.

See [docs/methods.md](docs/methods.md) and
[docs/input-formats.md](docs/input-formats.md) for details.

## Citation

Before publication, replace the placeholders in `CITATION.cff` with the author
list, article DOI, and final GitHub repository URL. Archive a release in Zenodo
to obtain a software DOI.

## License

MIT. See [LICENSE](LICENSE).
