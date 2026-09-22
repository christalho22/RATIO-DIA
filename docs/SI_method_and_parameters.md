# Supplementary Methods and Parameter Information for RATIO-DIA

## S1. Input data and initial spectral graph

Centroided DIA-MS/MS spectra were exported in MGF format, with the `SCANS`
identifier matched to the feature `Alignment ID` in the abundance table. After
feature alignment, 924 precursor features were retained. Pairwise spectral
similarity was evaluated using modified cosine matching with a fragment mass
tolerance of 0.01 Da. An original spectral edge was retained when the modified
cosine score was at least 0.70 and at least six fragment-ion pairs were matched.
This produced 15,417 original similarity edges. The 878 nodes incident to at
least one retained edge were used as the fixed input for module construction;
the other 46 aligned features were recorded as isolated and were not introduced
into the partition optimization.

## S2. Orthogonal fragmentation and topology evidence

For each spectrum, product ions below 50 Da or above the precursor mass plus
1 Da were removed. Peaks with relative intensity below 0.5% of the spectrum
maximum were discarded. Fragment masses were discretized into 0.02-Da bins.

Single-fragment evidence was represented by binary TF-IDF vectors after
retaining at most the 100 most intense filtered peaks per spectrum. Fragment
bins occurring in fewer than two spectra or in more than 60% of the input
spectra were excluded. Pairwise cosine similarity between the resulting vectors
defined the single-fragment similarity.

Fragment-pair evidence was constructed from unordered pairs of co-occurring
fragment bins after retaining at most the 40 most intense filtered peaks per
spectrum. Fragment-pair features occurring in fewer than two spectra or in more
than 20% of the input spectra were excluded. Cosine similarity between the
fragment-pair TF-IDF vectors defined the co-occurrence similarity. In the
reported dataset, 3,090 single-fragment features and 79,759 fragment-pair
features passed these filters.

Local topological support was calculated for each original edge using a weighted
closed-neighborhood edge-clustering value (ECV). ECV values were converted to
percentile ranks, denoted by Qij, to place topology evidence on a common scale.

## S3. Evidence integration and motif-neighborhood augmentation

The integrated fragmentation-pattern similarity was calculated as:

`Mij = (1 - wp) × Singleij + wp × Pairij`,

where wp is the fragment-pair contribution. The base weight of an original
modified-cosine edge was:

`Wbase_ij = Sij [f0 + (1 - f0)Mij] [t0 + (1 - t0)Qij]`,

where Sij is the modified-cosine score, f0 is the fragmentation-evidence floor,
and t0 is the topology-evidence floor.

For every node, its top-ranked integrated fragmentation relationships were
selected when Mij exceeded the motif minimum. The node-wise selections were
combined by union; reciprocal neighborhood membership was not required. The
final weight was:

`Wfinal_ij = I[(i,j) in E0]Wbase_ij + I[(i,j) in EM] beta Mij`,

where E0 is the original modified-cosine edge set, EM is the motif-neighborhood
edge set, and beta is the motif scale. Consequently, an existing spectral edge
within EM received an additive motif contribution, whereas a previously absent
motif edge was assigned beta Mij and was not required to meet the modified-
cosine or matched-fragment threshold.

## S4. Multiresolution module partition and parameter optimization

The evidence-weighted graph was partitioned using weighted Louvain community
detection. A fixed grid of 1,008 candidate settings was evaluated. Module-size
constraints were applied during candidate-partition selection and not within the
Louvain optimization itself. Candidate settings were excluded when either the
largest module or the module containing the target feature ions exceeded 100
nodes.

Before structural identification, a predefined group of target feature ions,
including precursor ions at *m/z* 571.2803, 518.2184, and 559.2807 together with
related feature ions, was used to optimize network partitioning. Eligible
partitions were ranked lexicographically by: (1) the number of target feature
ions co-clustered in one module; (2) the size of that module, favoring
compactness; (3) the proportion of target-ion pairs assigned to the same module;
(4) the number of nodes assigned to modules containing more than 10 nodes; and
(5) weighted modularity. These feature ions defined the parameter-optimization
objective; their structural identities and bioactivities were not used to
construct graph edges or determine module boundaries.

### Fixed search grid and selected setting

| Parameter | Search values | Selected value |
|---|---|---|
| Fragment-pair contribution, wp | 0.25, 0.50, 0.75, 0.90 | 0.25 |
| Motif neighborhood size, K | 10, 15, 20 | 15 |
| Minimum motif similarity | 0.10, 0.15, 0.20 | 0.10 |
| Fragmentation floor, f0 | 0.10, 0.25 | 0.10 |
| Topology floor, t0 | fixed | 0.25 |
| Motif scale, beta | 0.80, 1.20 | 0.80 |
| Louvain resolution | 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00 | 2.25 |
| Maximum eligible module size | fixed | 100 nodes |
| Random seed | fixed | 20260914 |

The selected setting generated 20,937 evidence-weighted relationships and 21
modules containing at least three nodes. All 878 connected input nodes were
assigned to these modules. Eighteen modules contained more than 10 nodes and
collectively included 855 nodes. The largest module contained 96 nodes, and 11
of the 12 target feature ions were co-clustered in M1.

## S5. Post hoc activity association

Oil Red O lipid accumulation was used as the fraction-level activity response.
For fraction g, the efficacy recovery index was:

`Eg = (mean Model - mean fraction g) / (mean Model - mean Control)`.

Feature abundances within each fraction were normalized to a total of 1,000,000
and transformed using `log1p` for trend analysis. Pearson, Spearman, and Kendall
associations with the five-fraction efficacy profile were exported as descriptive
measures. Candidate prioritization combined the target-fraction abundance
percentile, enrichment relative to the other fractions, target-fraction abundance
share, and agreement with the observed activity direction. Activity information
was introduced only after the modules had been fixed and therefore did not
affect graph construction or partitioning.

## S6. Software and reproducibility outputs

The workflow was implemented in Python 3.9 or later using NumPy 1.23 or later
and NetworkX 3.2 or later. The public implementation reports version 0.5.0.
Each run exports a JSON parameter record, the thresholded original edge table,
the final node and relation tables, module summaries, post hoc activity tables,
and the complete parameter-search table when parameter optimization is enabled. The final
relation table contains the modified-cosine score, matched-fragment count, ECV,
ECV percentile, single-fragment similarity, fragment-pair similarity, integrated
fragmentation similarity, original-edge base weight, additive motif weight,
final weight, edge origin, and module assignments.

## S7. Recommended reproducible command

```text
ratio-dia --mgf spectra.mgf --edges-in modified_cosine_edges.csv
  --abundance-table abundance.csv --activity-table activity.csv
  --partition-method relation-aware --input-scope edge-connected
  --min-cosine 0.70 --min-matched-peaks 6 --fragment-tolerance 0.01
  --relative-intensity-cutoff 0.005 --single-top-peaks 100
  --pair-top-peaks 40 --fragment-bin-width 0.02
  --pair-motif-weight 0.25 --motif-k 15 --motif-minimum 0.10
  --fragmentation-floor 0.10 --topology-floor 0.25
  --motif-scale 0.80 --louvain-resolution 2.25
  --min-module-size 3 --random-seed 20260914 --output output_directory
```

To repeat parameter optimization rather than applying the selected setting
directly, add `--parameter-search --target-feature-nodes target_feature_nodes.txt`.
The target-feature file should contain one Alignment ID per line.

## S8. Interpretation limits

The modules are evidence-integrated spectral communities rather than confirmed
compound families. Co-fragmentation, adducts, in-source fragments, non-covalent
aggregates, and co-eluting isomers may still produce redundant or ambiguous
nodes. Candidate-structure consistency and fraction-level activity association
support prioritization but do not establish compound identity or causal activity;
orthogonal structural and biological validation remains required.
