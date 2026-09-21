#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RATIO-DIA: relation-aware organization of DIA-MS/MS networks.

Pipeline
--------
MGF -> modified-cosine edges -> fragment/motif evidence integration
    -> multiresolution Louvain modules -> post hoc activity prioritization.

The validated default combines modified-cosine similarity, weighted local
edge-clustering support, single-fragment TF-IDF, and fragment-pair co-occurrence
TF-IDF. Motif-neighborhood relationships are defined from the union of node-wise
high-ranking fragmentation-pattern similarities; reciprocal membership is not
required. A fixed parameter grid can be calibrated with representative features,
but fraction-level activity labels are introduced only after module construction.
Results are *relation-aware spectral modules*, not confirmed compound families.

The implementation is self-contained apart from NumPy. Existing edge CSV files
can be reused with ``--edges-in``, which is recommended for parameter scans.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import networkx as nx

__version__ = "0.5.0"


@dataclass
class Spectrum:
    scan: int
    precursor_mz: float
    rt_min: float
    mz: np.ndarray
    intensity: np.ndarray


def _number(text, default=0.0):
    try:
        return float(str(text).split()[0])
    except (TypeError, ValueError, IndexError):
        return default


def read_mgf(path):
    """Read the fields needed by this workflow without loading extra packages."""
    spectra, params, peaks, in_block = [], {}, [], False

    def finish(index):
        if not peaks:
            return
        scan = int(_number(params.get("SCANS", index), index))
        precursor = _number(params.get("PEPMASS", 0.0))
        if "RTINMINUTES" in params:
            rt = _number(params["RTINMINUTES"])
        else:
            rt = _number(params.get("RTINSECONDS", 0.0)) / 60.0
        arr = np.asarray(peaks, dtype=float)
        good = np.isfinite(arr).all(axis=1) & (arr[:, 0] > 0) & (arr[:, 1] >= 0)
        arr = arr[good]
        if not len(arr):
            return
        arr = arr[np.argsort(arr[:, 0])]
        norm = np.linalg.norm(arr[:, 1])
        intensity = arr[:, 1] / norm if norm else arr[:, 1]
        spectra.append(Spectrum(scan, precursor, rt, arr[:, 0], intensity))

    with open(path, "r", encoding="utf-8-sig", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if line == "BEGIN IONS":
                params, peaks, in_block = {}, [], True
            elif line == "END IONS":
                finish(len(spectra) + 1)
                in_block = False
            elif in_block and "=" in line:
                key, value = line.split("=", 1)
                params[key.upper()] = value
            elif in_block and line:
                fields = line.split()
                if len(fields) >= 2:
                    try:
                        peaks.append((float(fields[0]), float(fields[1])))
                    except ValueError:
                        pass
    return spectra


def modified_cosine(a, b, tolerance):
    """Greedy one-to-one direct/precursor-shifted peak matching."""
    shift = b.precursor_mz - a.precursor_mz
    candidates = []
    for i, (mz_a, int_a) in enumerate(zip(a.mz, a.intensity)):
        for center in (mz_a, mz_a + shift):
            left = np.searchsorted(b.mz, center - tolerance, side="left")
            right = np.searchsorted(b.mz, center + tolerance, side="right")
            for j in range(left, right):
                candidates.append((int_a * b.intensity[j], i, j))
    candidates.sort(reverse=True)
    used_a, used_b, score, matches = set(), set(), 0.0, 0
    for product, i, j in candidates:
        if i not in used_a and j not in used_b:
            used_a.add(i)
            used_b.add(j)
            score += product
            matches += 1
    return min(1.0, max(0.0, score)), matches


EDGE_FIELDS = ["Source_Scan", "Target_Scan", "Source_Precursor_mz",
               "Target_Precursor_mz", "Precursor_Delta_mz", "Matched_Peaks",
               "Modified_Cosine_Score"]


def build_edges(spectra, threshold, min_matches, tolerance, output_path):
    edges = []
    with open(output_path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=EDGE_FIELDS)
        writer.writeheader()
        for i, source in enumerate(spectra[:-1]):
            for target in spectra[i + 1:]:
                score, matches = modified_cosine(source, target, tolerance)
                if score >= threshold and matches >= min_matches:
                    row = {
                        "Source_Scan": source.scan, "Target_Scan": target.scan,
                        "Source_Precursor_mz": source.precursor_mz,
                        "Target_Precursor_mz": target.precursor_mz,
                        "Precursor_Delta_mz": target.precursor_mz-source.precursor_mz,
                        "Matched_Peaks": matches,
                        "Modified_Cosine_Score": score,
                    }
                    writer.writerow(row)
                    edges.append((source.scan, target.scan, score, matches))
    return edges


def read_edges(path, threshold, min_matches):
    edges = []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            u = int(float(row.get("Source_Scan", row.get("source"))))
            v = int(float(row.get("Target_Scan", row.get("target"))))
            w = float(row.get("Modified_Cosine_Score",
                              row.get("score", row.get("similarity", 1))))
            m = int(float(row.get("Matched_Peaks", row.get("matched_peaks", 0))))
            if u != v and w >= threshold and m >= min_matches:
                edges.append((u, v, w, m))
    return edges


def top_k_filter(edges, k=0, mode="mutual"):
    """Sparsify a similarity graph using node-wise strongest neighbors.

    mutual: retain an edge only if each endpoint ranks the other in its top K;
    union:  retain it if either endpoint does.  Set K=0 to disable.
    """
    if not k or k < 1:
        return edges
    incident = defaultdict(list)
    for idx, (u, v, score, matches) in enumerate(edges):
        incident[u].append((score, matches, -v, idx))
        incident[v].append((score, matches, -u, idx))
    selected = {}
    for node, candidates in incident.items():
        selected[node] = {x[3] for x in sorted(candidates, reverse=True)[:k]}
    kept = []
    for idx, edge in enumerate(edges):
        u, v = edge[:2]
        decision = idx in selected[u] and idx in selected[v]
        if mode == "union":
            decision = idx in selected[u] or idx in selected[v]
        if decision:
            kept.append(edge)
    return kept


def graph_from_edges(nodes, edges):
    adjacency = {node: {} for node in nodes}
    matches = {}
    for u, v, weight, nmatch in edges:
        if v not in adjacency[u] or weight > adjacency[u][v]:
            adjacency[u][v] = adjacency[v][u] = weight
            matches[tuple(sorted((u, v)))] = nmatch
    return adjacency, matches


def edge_clustering_values(adjacency):
    """Weighted closed-neighborhood edge clustering values.

    Closed neighborhoods retain the direct u-v evidence even when an edge is
    not part of a triangle; this is stated explicitly because HC-PIN variants
    differ in their neighborhood convention.
    """
    strengths = {u: sum(nbrs.values()) + 1.0 for u, nbrs in adjacency.items()}
    scores = {}
    for u, nbrs in adjacency.items():
        for v, direct_w in nbrs.items():
            if u >= v:
                continue
            common = (set(adjacency[u]) | {u}) & (set(adjacency[v]) | {v})
            su = sv = 0.0
            for k in common:
                su += 1.0 if k == u else adjacency[u].get(k, direct_w if k == v else 0.0)
                sv += 1.0 if k == v else adjacency[v].get(k, direct_w if k == u else 0.0)
            scores[(u, v)] = (su * sv) / (strengths[u] * strengths[v])
    return scores


def module_condition(cluster, adjacency, mode, threshold):
    cluster = set(cluster)
    per_node = []
    for u in cluster:
        internal = sum(w for v, w in adjacency[u].items() if v in cluster)
        external = sum(w for v, w in adjacency[u].items() if v not in cluster)
        per_node.append((internal, external))
    if mode == "strong":
        return all(inside > threshold * outside for inside, outside in per_node)
    return sum(x for x, _ in per_node) > threshold * sum(y for _, y in per_node)


def hcpin_partition(adjacency, mode="weak", threshold=1.0, min_size=3,
                    selection_mode="quality-cut", resolution=1.0):
    """HC-PIN-style hierarchy, retaining lambda-modules along the merge trace.

    Edges are added in decreasing edge-clustering order.  Importantly, a merge
    is not rejected merely because an early two-node cluster is not yet a
    lambda-module.  Qualifying clusters are recorded as the hierarchy grows.
    The ECV-driven merge tree is cut by a weighted-modularity objective.  Each
    eligible parent lambda-module is compared with the best combination of its
    eligible descendants.  This avoids reducing the result to ordinary graph
    connected components, while still producing a non-overlapping Cytoscape
    module assignment.  ``largest-first`` is retained only to reproduce legacy
    outputs.
    """
    ecv = edge_clustering_values(adjacency)
    parent = {u: u for u in adjacency}
    members = {u: {u} for u in adjacency}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    candidates = []
    tree_members = {i: set(members[u]) for i, u in enumerate(adjacency)}
    tree_children = {}
    tree_eligible = {i: False for i in tree_members}
    root_tree = {u: i for i, u in enumerate(adjacency)}
    next_tree = len(tree_members)
    for (u, v), merge_ecv in sorted(ecv.items(), key=lambda item: (-item[1], item[0])):
        ru, rv = find(u), find(v)
        if ru == rv:
            continue
        left_tree, right_tree = root_tree[ru], root_tree[rv]
        candidate = members[ru] | members[rv]
        if len(members[ru]) < len(members[rv]):
            ru, rv = rv, ru
        parent[rv] = ru
        members[ru] = candidate
        del members[rv]
        del root_tree[rv]
        this_tree = next_tree
        next_tree += 1
        tree_members[this_tree] = set(candidate)
        tree_children[this_tree] = (left_tree, right_tree)
        eligible = (len(candidate) >= min_size and
                    module_condition(candidate, adjacency, mode, threshold))
        tree_eligible[this_tree] = eligible
        root_tree[ru] = this_tree
        if eligible:
            internal = [adjacency[x][y] for x in candidate for y in adjacency[x]
                        if x < y and y in candidate]
            candidates.append((this_tree, set(candidate), merge_ecv,
                               float(np.mean(internal)) if internal else 0.0))

    assignment = {u: 0 for u in adjacency}
    if selection_mode == "largest-first":
        candidates.sort(key=lambda x: (-len(x[1]), -x[2], -x[3], min(x[1])))
        selected, claimed = [], set()
        for _, cluster, _, _ in candidates:
            available = cluster - claimed
            if len(available) >= min_size:
                selected.append(set(available))
                claimed.update(available)
    else:
        total_weight = sum(w for u in adjacency for v, w in adjacency[u].items()
                           if u < v)
        degree_weight = {u: sum(adjacency[u].values()) for u in adjacency}

        def modularity_contribution(cluster):
            if total_weight <= 0:
                return float("-inf")
            internal = sum(adjacency[u][v] for u in cluster for v in adjacency[u]
                           if u < v and v in cluster)
            volume = sum(degree_weight[u] for u in cluster)
            return (internal/total_weight
                    - resolution*(volume/(2*total_weight))**2)

        memo = {}

        def best_cut(tree_id):
            if tree_id in memo:
                return memo[tree_id]
            child_score, child_modules = 0.0, []
            for child in tree_children.get(tree_id, ()):
                score, groups = best_cut(child)
                child_score += score
                child_modules.extend(groups)
            own_score = (modularity_contribution(tree_members[tree_id])
                         if tree_eligible.get(tree_id, False) else float("-inf"))
            # Prefer descendants on exact ties so a component root cannot
            # suppress an equally supported internal partition.
            if own_score > child_score + 1e-12:
                result = (own_score, [set(tree_members[tree_id])])
            else:
                result = (child_score, child_modules)
            memo[tree_id] = result
            return result

        selected = []
        for tree_id in root_tree.values():
            _, groups = best_cut(tree_id)
            selected.extend(groups)

    selected.sort(key=lambda cluster: (-len(cluster), min(cluster)))
    for module_id, cluster in enumerate(selected, 1):
        for node in cluster:
            assignment[node] = module_id
    return assignment, ecv


def ecv_consensus_partition(adjacency, min_size=3, ecv_power=2.0,
                            base_resolution=2.0, final_resolution=2.5,
                            consensus_runs=20, final_seed=20260905,
                            edge_order=None, consensus_seed_start=0):
    """Stable full-graph decomposition without deleting Top-K edges.

    Every threshold-qualified edge is retained. Its spectral weight is
    continuously modulated by the percentile rank of its ECV. Repeated Louvain
    partitions estimate an edge coassignment probability, which is used as a
    second continuous weight before the final multiresolution partition.
    Neither structural annotations nor activity labels enter the algorithm.
    """
    ecv = edge_clustering_values(adjacency)
    # Preserve the deterministic source/target order of the edge table when it
    # is available.  This exactly reproduces parameter-search runs when equal
    # ECV values occur; otherwise endpoint order is a deterministic fallback.
    edge_keys = ([tuple(sorted(key)) for key in edge_order]
                 if edge_order is not None else sorted(ecv))
    if len(edge_keys) != len(ecv) or set(edge_keys) != set(ecv):
        raise ValueError("edge_order must contain every undirected edge once")
    ecv_values = np.asarray([ecv[key] for key in edge_keys], dtype=float)
    ranks = np.empty(len(edge_keys), dtype=float)
    ranks[np.argsort(ecv_values, kind="mergesort")] = (
        np.arange(len(edge_keys)) + 1
    ) / max(1, len(edge_keys))

    base_graph = nx.Graph()
    base_graph.add_nodes_from(adjacency)
    base_weight = {}
    for index, (u, v) in enumerate(edge_keys):
        weight = adjacency[u][v] * (0.10 + ranks[index]) ** ecv_power
        base_weight[(u, v)] = weight
        base_graph.add_edge(u, v, weight=weight)

    def flat_assignment(groups):
        retained = [set(group) for group in groups if len(group) >= min_size]
        retained.sort(key=lambda group: (-len(group), min(group)))
        result = {node: 0 for node in adjacency}
        for module_id, group in enumerate(retained, 1):
            for node in group:
                result[node] = module_id
        return result

    runs = []
    for seed in range(consensus_seed_start,
                      consensus_seed_start + consensus_runs):
        groups = nx.community.louvain_communities(
            base_graph, weight="weight", resolution=base_resolution, seed=seed
        )
        runs.append(flat_assignment(groups))

    consensus_graph = nx.Graph()
    consensus_graph.add_nodes_from(adjacency)
    evidence = {}
    for index, (u, v) in enumerate(edge_keys):
        probability = sum(
            run[u] > 0 and run[u] == run[v] for run in runs
        ) / consensus_runs
        final_weight = base_weight[(u, v)] * (0.05 + probability) ** 2
        consensus_graph.add_edge(u, v, weight=final_weight)
        evidence[(u, v)] = {
            "ECV_Percentile": ranks[index],
            "Base_Reliability_Weight": base_weight[(u, v)],
            "Consensus_Coassignment": probability,
            "Final_Consensus_Weight": final_weight,
        }
    groups = nx.community.louvain_communities(
        consensus_graph, weight="weight", resolution=final_resolution,
        seed=final_seed
    )
    return flat_assignment(groups), ecv, evidence


def single_fragment_tfidf_similarity(
        spectra, node_ids, relative_intensity_cutoff=0.005, top_peaks=100,
        bin_width=0.02, minimum_document_frequency=2,
        maximum_document_fraction=0.60):
    """Cosine similarity of binary single-fragment TF-IDF vectors.

    Product ions are selected by relative intensity, limited to the strongest
    peaks, discretized by ``bin_width``, and filtered by document frequency.
    The returned matrix follows ``node_ids`` order.
    """
    node_ids = list(node_ids)
    node_index = {node: index for index, node in enumerate(node_ids)}
    feature_sets = [set() for _ in node_ids]
    for spectrum in spectra:
        if spectrum.scan not in node_index or not len(spectrum.mz):
            continue
        relative = spectrum.intensity / max(float(np.max(spectrum.intensity)), 1e-15)
        valid = ((relative >= relative_intensity_cutoff) &
                 (spectrum.mz >= 50.0) &
                 (spectrum.mz <= spectrum.precursor_mz + 1.0))
        selected = np.where(valid)[0]
        if len(selected) > top_peaks:
            selected = selected[np.argsort(relative[selected])[-top_peaks:]]
        feature_sets[node_index[spectrum.scan]] = set(
            np.rint(spectrum.mz[selected] / bin_width).astype(int).tolist()
        )

    frequency = Counter(feature for features in feature_sets for feature in features)
    keep = {
        feature for feature, count in frequency.items()
        if count >= minimum_document_frequency and
        count <= maximum_document_fraction * len(node_ids)
    }
    feature_sets = [features & keep for features in feature_sets]
    weight2 = {
        feature: (math.log((len(node_ids) + 1) / (frequency[feature] + 1)) + 1) ** 2
        for feature in keep
    }
    norms = [math.sqrt(sum(weight2[feature] for feature in features)) or 1.0
             for features in feature_sets]
    similarity = np.zeros((len(node_ids), len(node_ids)), dtype=np.float32)
    for i, left in enumerate(feature_sets):
        for j in range(i + 1, len(node_ids)):
            shared = left & feature_sets[j]
            if shared:
                score = sum(weight2[feature] for feature in shared) / (norms[i] * norms[j])
                similarity[i, j] = similarity[j, i] = score
    return similarity, len(keep)


def fragment_pair_tfidf_similarity(
        spectra, node_ids, relative_intensity_cutoff=0.005, top_peaks=40,
        bin_width=0.02, minimum_document_frequency=2,
        maximum_document_fraction=0.20):
    """Cosine similarity of TF-IDF weighted co-occurring fragment-bin pairs."""
    node_ids = list(node_ids)
    node_index = {node: index for index, node in enumerate(node_ids)}
    feature_sets = [set() for _ in node_ids]
    for spectrum in spectra:
        if spectrum.scan not in node_index or not len(spectrum.mz):
            continue
        relative = spectrum.intensity / max(float(np.max(spectrum.intensity)), 1e-15)
        valid = ((relative >= relative_intensity_cutoff) &
                 (spectrum.mz >= 50.0) &
                 (spectrum.mz <= spectrum.precursor_mz + 1.0))
        selected = np.where(valid)[0]
        if len(selected) > top_peaks:
            selected = selected[np.argsort(relative[selected])[-top_peaks:]]
        bins = sorted(set(np.rint(spectrum.mz[selected] / bin_width)
                          .astype(int).tolist()))
        feature_sets[node_index[spectrum.scan]] = set(combinations(bins, 2))

    postings = defaultdict(list)
    for row_index, features in enumerate(feature_sets):
        for feature in features:
            postings[feature].append(row_index)
    postings = {
        feature: rows for feature, rows in postings.items()
        if len(rows) >= minimum_document_frequency and
        len(rows) <= maximum_document_fraction * len(node_ids)
    }
    weight2 = {
        feature: (math.log((len(node_ids) + 1) / (len(rows) + 1)) + 1) ** 2
        for feature, rows in postings.items()
    }
    totals = np.zeros(len(node_ids), dtype=np.float64)
    overlap = np.zeros((len(node_ids), len(node_ids)), dtype=np.float32)
    for feature, rows in postings.items():
        value = weight2[feature]
        for i in rows:
            totals[i] += value
        for i, j in combinations(rows, 2):
            overlap[i, j] += value
            overlap[j, i] += value
    denominator = np.sqrt(totals[:, None] * totals[None, :])
    similarity = np.divide(
        overlap, denominator, out=np.zeros_like(overlap), where=denominator > 0
    )
    np.fill_diagonal(similarity, 0)
    return similarity, len(postings)


def _edge_percentiles(adjacency, edge_order):
    ecv = edge_clustering_values(adjacency)
    keys = [tuple(sorted(edge)) for edge in edge_order]
    if len(keys) != len(ecv) or set(keys) != set(ecv):
        raise ValueError("edge_order must contain every undirected original edge once")
    values = np.asarray([ecv[key] for key in keys], dtype=float)
    ranks = np.empty(len(keys), dtype=float)
    ranks[np.argsort(values, kind="mergesort")] = (
        np.arange(len(keys)) + 1
    ) / max(1, len(keys))
    return ecv, {key: float(ranks[index]) for index, key in enumerate(keys)}


def prepare_relation_evidence(
        node_ids, spectra, edges, relative_intensity_cutoff=0.005,
        single_top_peaks=100, pair_top_peaks=40, bin_width=0.02):
    """Calculate the invariant evidence layers once for a parameter scan."""
    node_ids = sorted(node_ids)
    node_index = {node: index for index, node in enumerate(node_ids)}
    edge_order = [(u, v) for u, v, _, _ in edges]
    adjacency, _ = graph_from_edges(node_ids, edges)
    ecv, ecv_percentile = _edge_percentiles(adjacency, edge_order)
    single, single_feature_count = single_fragment_tfidf_similarity(
        spectra, node_ids, relative_intensity_cutoff, single_top_peaks,
        bin_width=bin_width
    )
    pair, pair_feature_count = fragment_pair_tfidf_similarity(
        spectra, node_ids, relative_intensity_cutoff, pair_top_peaks,
        bin_width=bin_width
    )
    return {
        "node_ids": node_ids,
        "node_index": node_index,
        "ecv": ecv,
        "ecv_percentile": ecv_percentile,
        "single": single,
        "pair": pair,
        "Single_Fragment_Features": single_feature_count,
        "Fragment_Pair_Features": pair_feature_count,
    }


def build_relation_aware_graph(
        node_ids, spectra, edges, pair_motif_weight=0.25, motif_k=15,
        motif_minimum=0.10, fragmentation_floor=0.10,
        topology_floor=0.25, motif_scale=0.80,
        relative_intensity_cutoff=0.005, single_top_peaks=100,
        pair_top_peaks=40, bin_width=0.02, prepared=None):
    """Construct the evidence-weighted graph used in the reported analysis.

    Original modified-cosine edges are retained and continuously weighted by
    fragmentation-pattern and ECV evidence. Motif-neighborhood edges are the
    union of each node's strongest relationships and do not require reciprocal
    membership. Motif evidence is added to an existing edge or creates a new
    edge when no original modified-cosine edge exists.
    """
    if prepared is None:
        prepared = prepare_relation_evidence(
            node_ids, spectra, edges, relative_intensity_cutoff,
            single_top_peaks, pair_top_peaks, bin_width
        )
    node_ids = prepared["node_ids"]
    node_index = prepared["node_index"]
    ecv = prepared["ecv"]
    ecv_percentile = prepared["ecv_percentile"]
    single = prepared["single"]
    pair = prepared["pair"]
    motif = (1 - pair_motif_weight) * single + pair_motif_weight * pair

    graph = nx.Graph()
    graph.add_nodes_from(node_ids)
    evidence = {}
    for u, v, score, matched_peaks in edges:
        if u not in node_index or v not in node_index:
            continue
        i, j = node_index[u], node_index[v]
        key = tuple(sorted((u, v)))
        motif_score = float(motif[i, j])
        topology_score = ecv_percentile[key]
        base_weight = (
            score *
            (fragmentation_floor + (1 - fragmentation_floor) * motif_score) *
            (topology_floor + (1 - topology_floor) * topology_score)
        )
        attributes = {
            "weight": base_weight,
            "Modified_Cosine_Score": score,
            "Matched_Peaks": matched_peaks,
            "Edge_Clustering_Value": ecv[key],
            "ECV_Percentile": topology_score,
            "Single_Fragment_Similarity": float(single[i, j]),
            "Fragment_Pair_Similarity": float(pair[i, j]),
            "Integrated_Fragmentation_Similarity": motif_score,
            "Base_Original_Edge_Weight": base_weight,
            "Motif_Addition": 0.0,
            "Edge_Origin": "modified_cosine",
        }
        graph.add_edge(u, v, **attributes)
        evidence[key] = attributes

    motif_pairs = set()
    limit = min(max(0, motif_k), max(0, len(node_ids) - 1))
    if limit:
        for i in range(len(node_ids)):
            order = np.argsort(motif[i], kind="mergesort")[-limit:]
            for j in order:
                j = int(j)
                if i != j and motif[i, j] >= motif_minimum:
                    motif_pairs.add((min(i, j), max(i, j)))
    for i, j in sorted(motif_pairs):
        u, v = node_ids[i], node_ids[j]
        key = tuple(sorted((u, v)))
        addition = motif_scale * float(motif[i, j])
        if graph.has_edge(u, v):
            graph[u][v]["weight"] += addition
            graph[u][v]["Motif_Addition"] = addition
            graph[u][v]["Edge_Origin"] = "modified_plus_motif"
        else:
            attributes = {
                "weight": addition,
                "Modified_Cosine_Score": 0.0,
                "Matched_Peaks": 0,
                "Edge_Clustering_Value": "",
                "ECV_Percentile": 0.0,
                "Single_Fragment_Similarity": float(single[i, j]),
                "Fragment_Pair_Similarity": float(pair[i, j]),
                "Integrated_Fragmentation_Similarity": float(motif[i, j]),
                "Base_Original_Edge_Weight": 0.0,
                "Motif_Addition": addition,
                "Edge_Origin": "motif_neighborhood",
            }
            graph.add_edge(u, v, **attributes)
        evidence[key] = dict(graph[u][v])
    return graph, ecv, evidence, {
        "Single_Fragment_Features": prepared["Single_Fragment_Features"],
        "Fragment_Pair_Features": prepared["Fragment_Pair_Features"],
    }


def relation_partition(graph, resolution=2.25, seed=20260914,
                       min_module_size=3):
    """Weighted Louvain partition with deterministic module numbering."""
    groups = nx.community.louvain_communities(
        graph, weight="weight", resolution=resolution, seed=seed
    )
    groups = sorted((set(group) for group in groups),
                    key=lambda group: (-len(group), min(group)))
    assignment = {node: 0 for node in graph}
    module_id = 0
    for group in groups:
        if len(group) < min_module_size:
            continue
        module_id += 1
        for node in group:
            assignment[node] = module_id
    return groups, assignment


def partition_metrics(groups, node_ids, representative_nodes=()):
    """Summarize a candidate partition using the manuscript selection rules."""
    groups = sorted((set(group) for group in groups),
                    key=lambda group: (-len(group), min(group)))
    formal = [group for group in groups if len(group) >= 3]
    full_assignment = {
        node: index for index, group in enumerate(groups, 1) for node in group
    }
    representatives = set(representative_nodes) & set(node_ids)
    counts = Counter(full_assignment[node] for node in representatives)
    if counts:
        best_module, best_count = counts.most_common(1)[0]
        best_size = len(groups[best_module - 1])
        same_pairs = sum(value * (value - 1) // 2 for value in counts.values())
        possible_pairs = len(representatives) * (len(representatives) - 1) // 2
    else:
        best_module = best_count = best_size = same_pairs = possible_pairs = 0
    return {
        "Modules_ge3": len(formal),
        "Assigned_Nodes": sum(map(len, formal)),
        "Coverage": sum(map(len, formal)) / max(1, len(node_ids)),
        "Modules_gt10": sum(len(group) > 10 for group in formal),
        "Nodes_in_Modules_gt10": sum(len(group) for group in formal if len(group) > 10),
        "Largest_Module": len(groups[0]) if groups else 0,
        "Representative_Best_Module_Count": best_count,
        "Representative_Best_Module_Size": best_size,
        "Representative_Module_Count": len(counts),
        "Representative_Pair_Coclustering": (
            same_pairs / possible_pairs if possible_pairs else 0.0
        ),
        "Representative_Distribution": ";".join(
            map(str, sorted(counts.values(), reverse=True))
        ),
        "Representative_Best_Module_ID": best_module,
    }


def calibrate_relation_parameters(node_ids, spectra, edges, representatives,
                                  max_module_size=100, seed=20260914,
                                  relative_intensity_cutoff=0.005,
                                  single_top_peaks=100, pair_top_peaks=40,
                                  bin_width=0.02):
    """Run the fixed 1,008-candidate grid used for target-guided calibration."""
    representatives = set(representatives)
    if not representatives:
        raise ValueError("Parameter calibration requires representative nodes")
    prepared = prepare_relation_evidence(
        node_ids, spectra, edges, relative_intensity_cutoff,
        single_top_peaks, pair_top_peaks, bin_width
    )
    rows = []
    candidate_number = 0
    for pair_weight in (0.25, 0.50, 0.75, 0.90):
        for motif_k in (10, 15, 20):
            for motif_minimum in (0.10, 0.15, 0.20):
                for fragmentation_floor in (0.10, 0.25):
                    for motif_scale in (0.80, 1.20):
                        graph, _, _, feature_counts = build_relation_aware_graph(
                            node_ids, spectra, edges,
                            pair_motif_weight=pair_weight,
                            motif_k=motif_k,
                            motif_minimum=motif_minimum,
                            fragmentation_floor=fragmentation_floor,
                            topology_floor=0.25,
                            motif_scale=motif_scale,
                            prepared=prepared,
                        )
                        for resolution in (1.50, 1.75, 2.00, 2.25,
                                           2.50, 2.75, 3.00):
                            candidate_number += 1
                            groups, _ = relation_partition(
                                graph, resolution, seed, min_module_size=3
                            )
                            metrics = partition_metrics(groups, node_ids, representatives)
                            row = {
                                "Candidate": f"RATIO_{candidate_number:04d}",
                                "Pair_Motif_Weight": pair_weight,
                                "Motif_K": motif_k,
                                "Motif_Minimum": motif_minimum,
                                "Fragmentation_Floor": fragmentation_floor,
                                "Topology_Floor": 0.25,
                                "Motif_Scale": motif_scale,
                                "Resolution": resolution,
                                "Graph_Edges": graph.number_of_edges(),
                                "Weighted_Modularity": nx.community.modularity(
                                    graph, groups, weight="weight", resolution=resolution
                                ),
                                **feature_counts,
                                **metrics,
                            }
                            row["Eligible_Maximum_Size"] = int(
                                metrics["Largest_Module"] <= max_module_size and
                                metrics["Representative_Best_Module_Size"] <= max_module_size
                            )
                            rows.append(row)
    eligible = [row for row in rows if row["Eligible_Maximum_Size"]]
    if not eligible:
        raise RuntimeError("No candidate partition satisfied the module-size constraint")
    eligible.sort(key=lambda row: (
        -row["Representative_Best_Module_Count"],
        row["Representative_Best_Module_Size"],
        -row["Representative_Pair_Coclustering"],
        -row["Nodes_in_Modules_gt10"],
        -row["Weighted_Modularity"],
    ))
    rank = {row["Candidate"]: index for index, row in enumerate(eligible, 1)}
    for row in rows:
        row["Selection_Rank"] = rank.get(row["Candidate"], "")
    rows.sort(key=lambda row: (
        0 if row["Eligible_Maximum_Size"] else 1,
        int(row["Selection_Rank"]) if row["Selection_Rank"] else 10**9,
        row["Candidate"],
    ))
    return dict(eligible[0]), rows


def percentile_ranks(values):
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = (np.arange(len(values)) + 1) / max(1, len(values))
    return ranks


def _average_ranks(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return ranks


def _pearson(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def _spearman(x, y):
    return _pearson(_average_ranks(x), _average_ranks(y))


def _kendall_tau_b(x, y):
    concordant = discordant = tie_x = tie_y = 0
    for i in range(len(x)-1):
        for j in range(i+1, len(x)):
            dx, dy = np.sign(x[j]-x[i]), np.sign(y[j]-y[i])
            if dx == 0 and dy != 0: tie_x += 1
            elif dy == 0 and dx != 0: tie_y += 1
            elif dx*dy > 0: concordant += 1
            elif dx*dy < 0: discordant += 1
    den = math.sqrt((concordant+discordant+tie_x)*(concordant+discordant+tie_y))
    return (concordant-discordant)/den if den else np.nan


def activity_efficacy(path, group_column="Group", response_column="Response"):
    """Calculate recovery indices from a long-format biological response table."""
    grouped = defaultdict(list)
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.readline()
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        for row in csv.DictReader(handle, delimiter=delimiter):
            grouped[row[group_column].strip()].append(float(row[response_column]))
    required = {"Model", "Control", "A", "B", "C", "D", "E"}
    missing = required - set(grouped)
    if missing:
        raise ValueError(f"Activity table is missing groups: {sorted(missing)}")
    means = {g: float(np.mean(v)) for g, v in grouped.items()}
    den = means["Model"] - means["Control"]
    if den == 0:
        raise ValueError("Model and Control means are identical; recovery index is undefined")
    return {g: (means["Model"]-means[g])/den for g in "ABCDE"}


def read_abundance_and_label(path, efficacy):
    """Full A-E activity association and D-centric ranking."""
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.readline()
        delimiter = "\t" if "\t" in sample else ","
        handle.seek(0)
        rows = list(csv.DictReader(handle, delimiter=delimiter))
    parsed = []
    for row in rows:
        try:
            parsed.append({
                "id": int(float(row["Alignment ID"])),
                "rt": _number(row.get("Average Rt(min)", 0)),
                "mz": _number(row.get("Average Mz", 0)),
                **{g: _number(row.get(g, 0)) for g in "ABCDE"},
            })
        except (KeyError, ValueError):
            continue
    raw = np.asarray([[r[g] for g in "ABCDE"] for r in parsed], float)
    sums = raw.sum(axis=0)
    normalized = raw / sums * 1_000_000.0
    analysis = np.log1p(normalized)
    d_values = normalized[:, 3]
    pct = percentile_ranks(np.log1p(d_values))
    y = np.asarray([efficacy[g] for g in "ABCDE"])
    y_nod = np.asarray([efficacy[g] for g in "ABCE"])
    positives = normalized[normalized > 0]
    eps = max(1e-12, float(np.median(positives))*1e-6) if positives.size else 1e-12
    intermediate = []
    for row, norm, ana, p in zip(parsed, normalized, analysis, pct):
        abce_idx = [0, 1, 2, 4]
        abce = np.median(norm[abce_idx])
        log2fc = math.log2((norm[3] + eps) / (abce + eps))
        share = norm[3] / (norm.sum() + eps)
        pearson, spearman, kendall = _pearson(ana, y), _spearman(ana, y), _kendall_tau_b(ana, y)
        ana_nod = ana[abce_idx]
        pearson_nod, spearman_nod, kendall_nod = (_pearson(ana_nod, y_nod),
                                                  _spearman(ana_nod, y_nod),
                                                  _kendall_tau_b(ana_nod, y_nod))
        ordered = sorted("ABCDE", key=lambda g: efficacy[g])
        idx = {g: "ABCDE".index(g) for g in "ABCDE"}
        score_sum = weight_sum = 0.0
        for i, low in enumerate(ordered[:-1]):
            for high in ordered[i+1:]:
                weight = efficacy[high]-efficacy[low]
                delta = ana[idx[high]]-ana[idx[low]]
                score_sum += weight*(1.0 if delta > 1e-12 else (0.0 if delta < -1e-12 else .5))
                weight_sum += weight
        global_direction = score_sum/weight_sum
        dsum = dweight = 0.0
        for g in "ABCE":
            weight = efficacy["D"]-efficacy[g]
            delta = ana[3]-ana[idx[g]]
            dsum += weight*(1.0 if delta > 1e-12 else (0.0 if delta < -1e-12 else .5))
            dweight += weight
        d_direction = dsum/dweight
        direction = .65*d_direction + .35*global_direction
        intermediate.append({**row, **{f"{g}_Normalized": norm[i] for i,g in enumerate("ABCDE")},
                             "D_percentile": p, "D_log2FC_vs_ABCE": log2fc,
                             "D_share_AE": share, "D_direction": d_direction,
                             "Global_activity_direction": global_direction,
                             "Activity_direction_score": direction,
                             "Pearson_r_AE": pearson, "Spearman_rho_AE": spearman,
                             "Kendall_tau_AE": kendall, "Pearson_r_NoD": pearson_nod,
                             "Spearman_rho_NoD": spearman_nod, "Kendall_tau_NoD": kendall_nod})
    enrich_pct = percentile_ranks(np.asarray([r["D_log2FC_vs_ABCE"] for r in intermediate]))
    share_pct = percentile_ranks(np.asarray([r["D_share_AE"] for r in intermediate]))
    result = {}
    for row, ep, sp in zip(intermediate, enrich_pct, share_pct):
        specificity = .70*ep + .30*sp
        # Same fallback quality score as the original when GNPS table has no S/N, Fill%, MS/MS flag.
        quality = .35
        priority = 100*(.30*row["D_percentile"] + .35*specificity
                        + .25*row["Activity_direction_score"] + .10*quality)
        tier1 = (row["D_percentile"] >= .95 and row["D_log2FC_vs_ABCE"] >= 1
                 and row["D_share_AE"] >= .50 and row["Activity_direction_score"] >= .85)
        tier2 = (row["D_percentile"] >= .90 and row["D_log2FC_vs_ABCE"] >= 1
                 and row["D_share_AE"] >= .40 and not tier1)
        tier = 1 if tier1 else (2 if tier2 else 0)
        result[row["id"]] = {**row, "D_specificity_score": specificity,
                             "D_centric_priority_score": priority,
                             "D_tier": tier, "D_candidate": int(tier in (1, 2)),
                             **{f"Efficacy_{g}": efficacy[g] for g in "ABCDE"}}
    return result


def write_outputs(outdir, spectra, edges, adjacency, matches, assignment, ecv,
                  abundance, edge_evidence=None, weighted_graph=None):
    outdir.mkdir(parents=True, exist_ok=True)
    spec_by_id = {s.scan: s for s in spectra}
    activity_fields = ["D_candidate", "D_tier", "A", "B", "C", "D", "E", "D_percentile",
                       "D_log2FC_vs_ABCE", "D_share_AE", "D_direction",
                       "Global_activity_direction", "Activity_direction_score",
                       "Pearson_r_AE", "Spearman_rho_AE", "Kendall_tau_AE",
                       "Pearson_r_NoD", "Spearman_rho_NoD", "Kendall_tau_NoD",
                       "D_specificity_score", "D_centric_priority_score",
                       "Efficacy_A", "Efficacy_B", "Efficacy_C", "Efficacy_D", "Efficacy_E"]
    node_fields = ["Alignment_ID", "Precursor_mz", "RT_min", "Module_ID"] + activity_fields
    with open(outdir/"02_nodes_with_modules_and_D_labels.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=node_fields); w.writeheader()
        for node in sorted(adjacency):
            s, a = spec_by_id.get(node), abundance.get(node, {})
            w.writerow({"Alignment_ID": node, "Precursor_mz": s.precursor_mz if s else a.get("mz", ""),
                        "RT_min": s.rt_min if s else a.get("rt", ""), "Module_ID": assignment[node],
                        **{k: a.get(k, "") for k in activity_fields}})

    relation_fields = ["ECV_Percentile", "Single_Fragment_Similarity",
                       "Fragment_Pair_Similarity",
                       "Integrated_Fragmentation_Similarity",
                       "Base_Original_Edge_Weight", "Motif_Addition",
                       "Final_Weight", "Edge_Origin"]
    consensus_fields = ["ECV_Percentile", "Base_Reliability_Weight",
                        "Consensus_Coassignment", "Final_Consensus_Weight"]
    edge_fields = EDGE_FIELDS + ["Edge_Clustering_Value"]
    if weighted_graph is not None:
        edge_fields += relation_fields
    elif edge_evidence is not None:
        edge_fields += consensus_fields
    edge_fields += ["Source_Module", "Target_Module", "Same_Module"]
    with open(outdir/"03_edges_with_modules.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=edge_fields); w.writeheader()
        if weighted_graph is not None:
            edge_iterator = ((u, v, data) for u, v, data in weighted_graph.edges(data=True))
            for u, v, data in edge_iterator:
                su, sv = spec_by_id.get(u), spec_by_id.get(v)
                w.writerow({
                    "Source_Scan": u, "Target_Scan": v,
                    "Source_Precursor_mz": su.precursor_mz if su else "",
                    "Target_Precursor_mz": sv.precursor_mz if sv else "",
                    "Precursor_Delta_mz": (sv.precursor_mz-su.precursor_mz) if su and sv else "",
                    "Matched_Peaks": data.get("Matched_Peaks", 0),
                    "Modified_Cosine_Score": data.get("Modified_Cosine_Score", 0),
                    "Edge_Clustering_Value": data.get("Edge_Clustering_Value", ""),
                    "ECV_Percentile": data.get("ECV_Percentile", 0),
                    "Single_Fragment_Similarity": data.get("Single_Fragment_Similarity", 0),
                    "Fragment_Pair_Similarity": data.get("Fragment_Pair_Similarity", 0),
                    "Integrated_Fragmentation_Similarity": data.get("Integrated_Fragmentation_Similarity", 0),
                    "Base_Original_Edge_Weight": data.get("Base_Original_Edge_Weight", 0),
                    "Motif_Addition": data.get("Motif_Addition", 0),
                    "Final_Weight": data.get("weight", 0),
                    "Edge_Origin": data.get("Edge_Origin", ""),
                    "Source_Module": assignment[u], "Target_Module": assignment[v],
                    "Same_Module": int(assignment[u] > 0 and assignment[u] == assignment[v]),
                })
        else:
            for u, v, weight, nmatch in edges:
                su, sv = spec_by_id.get(u), spec_by_id.get(v)
                extra = ((edge_evidence or {}).get(tuple(sorted((u, v))), {}))
                w.writerow({"Source_Scan": u, "Target_Scan": v,
                            "Source_Precursor_mz": su.precursor_mz if su else "",
                            "Target_Precursor_mz": sv.precursor_mz if sv else "",
                            "Precursor_Delta_mz": (sv.precursor_mz-su.precursor_mz) if su and sv else "",
                            "Matched_Peaks": nmatch, "Modified_Cosine_Score": weight,
                            "Edge_Clustering_Value": ecv.get(tuple(sorted((u, v))), ""),
                            **extra,
                            "Source_Module": assignment[u], "Target_Module": assignment[v],
                            "Same_Module": int(assignment[u] > 0 and assignment[u] == assignment[v])})

    modules = defaultdict(list)
    for node, mid in assignment.items():
        modules[mid].append(node)
    summary_fields = ["Module_ID", "Node_Count", "Internal_Edge_Count", "Density",
                      "Mean_Internal_Cosine", "D_Candidate_Count", "D_Candidate_Ratio",
                      "Contains_D_Candidate", "Sum_A", "Sum_B", "Sum_C", "Sum_D", "Sum_E"]
    summaries = []
    for mid, nodes in sorted(modules.items()):
        node_set, internal = set(nodes), []
        for u in nodes:
            if weighted_graph is not None:
                internal += [data.get("Modified_Cosine_Score", 0)
                             for v, data in weighted_graph[u].items()
                             if u < v and v in node_set and data.get("Modified_Cosine_Score", 0) > 0]
            else:
                internal += [w for v, w in adjacency[u].items() if u < v and v in node_set]
        n, dcount = len(nodes), sum(abundance.get(x, {}).get("D_candidate", 0) for x in nodes)
        row = {"Module_ID": mid, "Node_Count": n, "Internal_Edge_Count": len(internal),
               "Density": 2*len(internal)/(n*(n-1)) if n > 1 else 0,
               "Mean_Internal_Cosine": np.mean(internal) if internal else 0,
               "D_Candidate_Count": dcount, "D_Candidate_Ratio": dcount/n if n else 0,
               "Contains_D_Candidate": int(dcount > 0)}
        row.update({f"Sum_{g}": sum(abundance.get(x, {}).get(g, 0) for x in nodes) for g in "ABCDE"})
        summaries.append(row)
    with open(outdir/"04_module_summary.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=summary_fields); w.writeheader(); w.writerows(summaries)

    focus_modules = {r["Module_ID"] for r in summaries if r["Module_ID"] > 0 and r["Contains_D_Candidate"]}
    focus_nodes = {n for n, mid in assignment.items() if mid in focus_modules}
    with open(outdir/"05_D_focused_nodes.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=node_fields); w.writeheader()
        for node in sorted(focus_nodes):
            s, a = spec_by_id.get(node), abundance.get(node, {})
            w.writerow({"Alignment_ID": node, "Precursor_mz": s.precursor_mz if s else "",
                        "RT_min": s.rt_min if s else "", "Module_ID": assignment[node],
                        **{k: a.get(k, "") for k in activity_fields}})
    with open(outdir/"06_D_focused_edges.csv", "w", newline="", encoding="utf-8-sig") as f:
        fields = ["Source_Scan", "Target_Scan", "Modified_Cosine_Score", "Matched_Peaks", "Module_ID"]
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        if weighted_graph is not None:
            focus_iterator = ((u, v, data.get("Modified_Cosine_Score", 0),
                               data.get("Matched_Peaks", 0))
                              for u, v, data in weighted_graph.edges(data=True))
        else:
            focus_iterator = edges
        for u, v, score, nmatch in focus_iterator:
            if u in focus_nodes and v in focus_nodes and assignment[u] == assignment[v]:
                w.writerow({"Source_Scan": u, "Target_Scan": v, "Modified_Cosine_Score": score,
                            "Matched_Peaks": nmatch, "Module_ID": assignment[u]})

    candidate_fields = ["Module_ID", "Alignment_ID", "Precursor_mz", "RT_min"] + activity_fields[1:] + ["Module_Node_Count",
                        "Module_D_Candidate_Count"]
    summary_by_module = {r["Module_ID"]: r for r in summaries}
    candidate_rows = []
    for node, a in abundance.items():
        mid = assignment.get(node, 0)
        if mid > 0 and a.get("D_candidate"):
            s = spec_by_id.get(node)
            sm = summary_by_module[mid]
            candidate_rows.append({"Module_ID": mid, "Alignment_ID": node,
                                   "Precursor_mz": s.precursor_mz if s else a.get("mz", ""),
                                   "RT_min": s.rt_min if s else a.get("rt", ""),
                                   **{k: a.get(k, "") for k in activity_fields[1:]},
                                   "Module_Node_Count": sm["Node_Count"],
                                   "Module_D_Candidate_Count": sm["D_Candidate_Count"]})
    candidate_rows.sort(key=lambda r: (r["Module_ID"], r["D_tier"],
                                       -r["D_centric_priority_score"]))
    with open(outdir/"07_D_active_marker_mz_by_module.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=candidate_fields); w.writeheader(); w.writerows(candidate_rows)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mgf", required=True)
    p.add_argument("--abundance-table", required=True)
    p.add_argument("--activity-table", required=True,
                   help="Long CSV/TSV with Group and Response columns; groups: Model, Control, A-E")
    p.add_argument("--activity-group-column", default="Group")
    p.add_argument("--activity-response-column", default="Response")
    p.add_argument("--edges-in", help="Reuse an existing edge CSV instead of recalculating spectra")
    p.add_argument("--output", default="ratio_dia_output")
    p.add_argument("--min-cosine", type=float, default=.70)
    p.add_argument("--min-matched-peaks", type=int, default=6)
    p.add_argument("--top-k", type=int, default=0,
                   help="Optional strongest-K prefilter; 0 (default) sends the complete thresholded graph directly to hierarchical partitioning")
    p.add_argument("--top-k-mode", choices=["mutual", "union"], default="mutual")
    p.add_argument("--fragment-tolerance", type=float, default=.01)
    p.add_argument("--mode", choices=["weak", "strong"], default="weak")
    p.add_argument("--hcpin-threshold", type=float, default=1.0)
    p.add_argument("--min-module-size", type=int, default=3)
    p.add_argument("--selection-mode", choices=["quality-cut", "largest-first"],
                   default="quality-cut",
                   help="Flat cut of the HC-PIN hierarchy; largest-first reproduces legacy output")
    p.add_argument("--resolution", type=float, default=1.0,
                   help="Weighted-modularity resolution for quality-cut")
    p.add_argument("--partition-method",
                   choices=["relation-aware", "ecv-consensus", "ecv-hierarchy"],
                   default="relation-aware",
                   help="relation-aware reproduces the reported multi-evidence workflow")
    p.add_argument("--input-scope", choices=["edge-connected", "all"],
                   default="edge-connected",
                   help="Use only nodes retained by the thresholded similarity graph (reported default), or all MGF spectra")
    p.add_argument("--pair-motif-weight", type=float, default=.25)
    p.add_argument("--motif-k", type=int, default=15)
    p.add_argument("--motif-minimum", type=float, default=.10)
    p.add_argument("--fragmentation-floor", type=float, default=.10)
    p.add_argument("--topology-floor", type=float, default=.25)
    p.add_argument("--motif-scale", type=float, default=.80)
    p.add_argument("--relative-intensity-cutoff", type=float, default=.005)
    p.add_argument("--single-top-peaks", type=int, default=100)
    p.add_argument("--pair-top-peaks", type=int, default=40)
    p.add_argument("--fragment-bin-width", type=float, default=.02)
    p.add_argument("--louvain-resolution", type=float, default=2.25)
    p.add_argument("--max-module-size", type=int, default=100)
    p.add_argument("--representative-nodes",
                   help="Comma-separated IDs or a one-column text file used only for calibration")
    p.add_argument("--parameter-search", action="store_true",
                   help="Run the fixed 1,008-candidate calibration grid")
    p.add_argument("--ecv-power", type=float, default=2.0)
    p.add_argument("--consensus-base-resolution", type=float, default=2.0)
    p.add_argument("--consensus-final-resolution", type=float, default=2.5)
    p.add_argument("--consensus-runs", type=int, default=20)
    p.add_argument("--random-seed", type=int, default=20260914)
    return p.parse_args()


def parse_representative_nodes(value):
    if not value:
        return []
    path = Path(value)
    if path.exists():
        text = path.read_text(encoding="utf-8-sig")
    else:
        text = value
    tokens = text.replace(",", " ").replace(";", " ").split()
    return [int(float(token)) for token in tokens]


def main():
    args = parse_args()
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Reading MGF: {args.mgf}")
    spectra = read_mgf(args.mgf)
    print(f"Spectra: {len(spectra)}")
    raw_edge_path = outdir/"01_filtered_similarity_edges.csv"
    if args.edges_in:
        edges = read_edges(args.edges_in, args.min_cosine, args.min_matched_peaks)
        with open(raw_edge_path, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f); w.writerow(EDGE_FIELDS)
            by_id = {s.scan: s for s in spectra}
            for u, v, score, nmatch in edges:
                su, sv = by_id.get(u), by_id.get(v)
                w.writerow([u, v, su.precursor_mz if su else "", sv.precursor_mz if sv else "",
                            sv.precursor_mz-su.precursor_mz if su and sv else "", nmatch, score])
    else:
        print("Calculating all pairwise modified-cosine similarities...")
        edges = build_edges(spectra, args.min_cosine, args.min_matched_peaks,
                            args.fragment_tolerance, raw_edge_path)
    pre_topk_edges = len(edges)
    edges = top_k_filter(edges, args.top_k, args.top_k_mode)
    if args.top_k:
        print(f"Top-K ({args.top_k_mode}, K={args.top_k}): {pre_topk_edges} -> {len(edges)} edges")
    connected_nodes = {node for edge in edges for node in edge[:2]}
    nodes = connected_nodes if args.input_scope == "edge-connected" else {s.scan for s in spectra}
    edges = [edge for edge in edges if edge[0] in nodes and edge[1] in nodes]
    selected_spectra = [spectrum for spectrum in spectra if spectrum.scan in nodes]
    adjacency, matches = graph_from_edges(nodes, edges)
    print(f"Input scope: {len(nodes)} nodes; retained edges: {len(edges)}")
    edge_evidence = None
    weighted_graph = None
    selected_parameters = {
        "pair_motif_weight": args.pair_motif_weight,
        "motif_k": args.motif_k,
        "motif_minimum": args.motif_minimum,
        "fragmentation_floor": args.fragmentation_floor,
        "topology_floor": args.topology_floor,
        "motif_scale": args.motif_scale,
        "resolution": args.louvain_resolution,
    }
    feature_counts = {}
    if args.partition_method == "relation-aware":
        if args.parameter_search:
            representatives = parse_representative_nodes(args.representative_nodes)
            best, search_rows = calibrate_relation_parameters(
                nodes, selected_spectra, edges, representatives,
                args.max_module_size, args.random_seed,
                args.relative_intensity_cutoff, args.single_top_peaks,
                args.pair_top_peaks, args.fragment_bin_width
            )
            with open(outdir/"08_parameter_search.csv", "w", newline="",
                      encoding="utf-8-sig") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(search_rows[0]))
                writer.writeheader(); writer.writerows(search_rows)
            selected_parameters = {
                "pair_motif_weight": best["Pair_Motif_Weight"],
                "motif_k": best["Motif_K"],
                "motif_minimum": best["Motif_Minimum"],
                "fragmentation_floor": best["Fragmentation_Floor"],
                "topology_floor": best["Topology_Floor"],
                "motif_scale": best["Motif_Scale"],
                "resolution": best["Resolution"],
            }
        weighted_graph, ecv, edge_evidence, feature_counts = build_relation_aware_graph(
            nodes, selected_spectra, edges,
            pair_motif_weight=selected_parameters["pair_motif_weight"],
            motif_k=selected_parameters["motif_k"],
            motif_minimum=selected_parameters["motif_minimum"],
            fragmentation_floor=selected_parameters["fragmentation_floor"],
            topology_floor=selected_parameters["topology_floor"],
            motif_scale=selected_parameters["motif_scale"],
            relative_intensity_cutoff=args.relative_intensity_cutoff,
            single_top_peaks=args.single_top_peaks,
            pair_top_peaks=args.pair_top_peaks,
            bin_width=args.fragment_bin_width,
        )
        groups, assignment = relation_partition(
            weighted_graph, selected_parameters["resolution"],
            args.random_seed, args.min_module_size
        )
    elif args.partition_method == "ecv-consensus":
        assignment, ecv, edge_evidence = ecv_consensus_partition(
            adjacency, args.min_module_size, args.ecv_power,
            args.consensus_base_resolution, args.consensus_final_resolution,
            args.consensus_runs, args.random_seed,
            [(u, v) for u, v, _, _ in edges])
    else:
        assignment, ecv = hcpin_partition(
            adjacency, args.mode, args.hcpin_threshold,
            args.min_module_size, args.selection_mode, args.resolution)
    efficacy = activity_efficacy(args.activity_table, args.activity_group_column,
                                 args.activity_response_column)
    print("Recovery indices:", ", ".join(f"{k}={v:.4f}" for k, v in efficacy.items()))
    abundance = read_abundance_and_label(args.abundance_table, efficacy)
    write_outputs(outdir, selected_spectra, edges, adjacency, matches, assignment, ecv,
                  abundance, edge_evidence, weighted_graph)
    run_parameters = {
        "software_version": __version__,
        "partition_method": args.partition_method,
        "input_scope": args.input_scope,
        "input_nodes": len(nodes),
        "original_similarity_edges": len(edges),
        "final_weighted_edges": weighted_graph.number_of_edges() if weighted_graph is not None else len(edges),
        "minimum_modified_cosine": args.min_cosine,
        "minimum_matched_fragments": args.min_matched_peaks,
        "fragment_tolerance_Da": args.fragment_tolerance,
        "relative_intensity_cutoff": args.relative_intensity_cutoff,
        "single_fragment_top_peaks": args.single_top_peaks,
        "fragment_pair_top_peaks": args.pair_top_peaks,
        "fragment_bin_width_Da": args.fragment_bin_width,
        "random_seed": args.random_seed,
        "parameter_search": args.parameter_search,
        "selected_parameters": selected_parameters,
        "evidence_feature_counts": feature_counts,
    }
    (outdir/"00_run_parameters.json").write_text(
        json.dumps(run_parameters, indent=2), encoding="utf-8"
    )
    module_count = len({x for x in assignment.values() if x > 0})
    focused = len({assignment[n] for n in assignment if assignment[n] > 0 and abundance.get(n, {}).get("D_candidate")})
    print(f"RATIO-DIA modules: {module_count}; modules containing D candidates: {focused}")
    print(f"Output: {outdir.resolve()}")


if __name__ == "__main__":
    main()
