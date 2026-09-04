# Method summary

## Spectral graph construction

Each aligned precursor feature is represented as a node. An undirected edge is
created when two spectra meet the modified-cosine threshold and minimum matched-
peak requirement. For spectra \(i\) and \(j\), normalized fragment intensities
are matched either directly or after applying their precursor-mass shift. The
edge weight is the sum of products of greedily selected one-to-one peak matches.

## Reciprocal-neighborhood focusing

For each node, incident edges are ranked by spectral similarity and matched-peak
count. An edge is retained in mutual mode only when both endpoints place one
another among their strongest \(K\) neighbors. This limits hub-driven and weak
cross-cluster bridging while preserving locally representative relationships.

## Topology-focused modules

Retained edges are ordered by a weighted closed-neighborhood edge-clustering
value. Connected sets are accumulated hierarchically and recorded when their
weighted internal connectivity exceeds \(\lambda\) times their external
connectivity. In weak mode this condition is evaluated at module level; in strong
mode it must hold for every node. The final export is a non-overlapping flat
partition. These are spectral-topology modules, not confirmed chemical families.

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

For every annotated node pair within a module, all CHON candidates are enumerated.
The workflow records the maximum Morgan-fingerprint Tanimoto similarity and
whether any candidate combination shares a Murcko scaffold or chemical ontology.
This avoids relying exclusively on Top-1 in-silico annotations, but it is a
permissive measure and must not be interpreted as structural confirmation.
