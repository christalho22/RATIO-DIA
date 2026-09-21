# Method summary

## Spectral graph construction

Each aligned precursor feature is represented as a node. An undirected edge is
created when two spectra meet the modified-cosine threshold and minimum matched-
peak requirement. For spectra \(i\) and \(j\), normalized fragment intensities
are matched either directly or after applying their precursor-mass shift. The
edge weight is the sum of products of greedily selected one-to-one peak matches.

## Optional reciprocal-neighborhood sensitivity analysis

The primary workflow applies no neighborhood prefilter: all edges satisfying the
spectral-similarity and matched-peak criteria are passed directly to the
hierarchical topology analysis. For sensitivity analysis only, incident edges
may be ranked by spectral similarity and matched-peak count and retained in
mutual mode when both endpoints place one another among their strongest \(K\)
neighbors. Results using this optional prefilter must not be presented as the
primary no-Top-K decomposition.

## Relation-aware evidence integration (primary analysis)

Three evidence layers are derived independently of biological activity. First,
the weighted closed-neighborhood edge-clustering value (ECV) of each original
modified-cosine edge is transformed to a percentile, \(Q_{ij}\). Second,
filtered product-ion bins are represented by binary TF-IDF vectors and compared
by cosine similarity. Third, unordered pairs of co-occurring product-ion bins
are represented by TF-IDF vectors. The single-fragment and fragment-pair
similarities are combined as \(M_{ij}\).

For an original spectral edge, the base evidence weight is

\[
W^{base}_{ij}=S_{ij}[f_0+(1-f_0)M_{ij}]
              [t_0+(1-t_0)Q_{ij}],
\]

where \(S_{ij}\) is the modified-cosine score and the positive floors retain
weak but non-zero orthogonal evidence. A motif neighborhood is formed from the
union of each node's top-ranked \(M_{ij}\) relationships above a minimum motif
similarity; reciprocal membership is not required. The final relation weight is

\[
W^{final}_{ij}=\mathbf{1}_{(i,j)\in E_0}W^{base}_{ij}
               +\mathbf{1}_{(i,j)\in E_M}\beta M_{ij}.
\]

Thus an original edge that also belongs to the motif neighborhood receives an
additive motif contribution, whereas a previously absent motif edge is assigned
\(\beta M_{ij}\) without being required to satisfy the modified-cosine cutoff.

## Multiresolution topology partition and calibration

Weighted Louvain partitions are generated over a fixed parameter grid. The
module-size constraint is applied during candidate-partition selection rather
than within Louvain optimization. Eligible partitions are ranked
lexicographically by: (1) the number of representative features recovered in a
single module, (2) the compactness of that target module, (3) representative-
pair co-clustering, (4) coverage of nodes in modules containing more than 10
nodes, and (5) weighted modularity. Representative features therefore serve as
a calibration set and not an independent validation set. Candidate structures,
fraction abundances, and activity measurements do not enter graph construction
or partition optimization.

For the reported dataset, the selected configuration used pair-motif weight
0.25, motif top-15 union neighborhoods, minimum motif similarity 0.10,
fragmentation and topology floors 0.10 and 0.25, motif scale 0.80, Louvain
resolution 2.25, and random seed 20260914. Exact preprocessing settings and the
full search grid are listed in the supplementary parameter document.

## Legacy ECV-consensus method (sensitivity analysis)

The earlier implementation, retained for reproducibility, continuously weights
original edges by ECV percentile, estimates coassignment across repeated
Louvain runs, and partitions the resulting consensus graph. It is not the
primary method used for the current manuscript.

## HC-PIN-inspired hierarchy (sensitivity analysis)

Retained edges are ordered by a weighted closed-neighborhood edge-clustering
value. Connected sets are accumulated hierarchically and recorded when their
weighted internal connectivity exceeds \(\lambda\) times their external
connectivity. In weak mode this condition is evaluated at module level; in strong
mode it must hold for every node. To obtain a non-overlapping flat partition,
each eligible parent is compared with the best combination of eligible descendants
using its weighted-modularity contribution,

\[
q(C)=\frac{w_{in}(C)}{W}-\gamma\left(\frac{s(C)}{2W}\right)^2,
\]

where \(w_{in}(C)\) is internal edge weight, \(s(C)\) is module strength,
\(W\) is total graph edge weight, and \(\gamma\) is the resolution parameter.
Dynamic programming selects the higher-scoring parent or descendant cut at each
branch. Consequently, the output does not collapse to connected components even
when an entire component satisfies the weak lambda condition. This quality cut
is a RATIO-DIA adaptation of the HC-PIN-style hierarchy, not the unmodified
original HC-PIN algorithm. It is retained as an alternative partition method.
Both outputs are spectral-topology modules, not confirmed chemical families.

## Activity association

Oil Red O response is summarized for Model, Control, and fractions A-E. The
recovery index of fraction \(g\) is

\[
E_g = \frac{\bar{Y}_{Model}-\bar{Y}_g}
           {\bar{Y}_{Model}-\bar{Y}_{Control}}.
\]

Feature abundances are normalized to a constant column sum and transformed with
`log1p`. Candidate prioritization combines target-fraction abundance percentile,
enrichment relative to the median of the other fractions, target-fraction share,
and directional agreement with the observed activity ordering. Pearson,
Spearman, and Kendall associations are exported as supporting descriptors but
are not treated as inferential tests because only five fraction-level points are
available.

The included Tier 1 and Tier 2 thresholds reproduce the study implementation.
They should be prespecified or sensitivity-tested before application to another
dataset.

## Structural consistency

For every annotated node pair within a module, all available candidates are enumerated.
The workflow records the maximum Morgan-fingerprint Tanimoto similarity and
whether any candidate combination shares a Murcko scaffold or chemical ontology.
This avoids relying exclusively on Top-1 in-silico annotations, but it is a
permissive measure and must not be interpreted as structural confirmation.
