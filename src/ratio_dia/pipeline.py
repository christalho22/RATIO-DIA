#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RATIO-DIA: activity-guided topology focusing of DIA-MS/MS networks.

Pipeline
--------
MGF -> modified-cosine edges -> HC-PIN-style agglomerative modules
    -> D-centric labels from the GNPS A-E abundance table -> focused D modules.

The HC-PIN implementation is an adaptation for an MS/MS similarity graph.  It
uses a weighted edge-clustering value to order candidate merges and accepts a
merge only when the resulting cluster satisfies a weak or strong lambda-module
condition.  Results should be described as *topological spectral modules*, not
as confirmed compounds or confirmed active molecular families.

The implementation is self-contained apart from NumPy. Existing edge CSV files
can be reused with ``--edges-in``, which is recommended for parameter scans.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__version__ = "0.1.0"


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
            w = float(row.get("Modified_Cosine_Score", row.get("score", 1)))
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


def hcpin_partition(adjacency, mode="weak", threshold=1.0, min_size=3):
    """HC-PIN-style hierarchy, retaining lambda-modules along the merge trace.

    Edges are added in decreasing edge-clustering order.  Importantly, a merge
    is not rejected merely because an early two-node cluster is not yet a
    lambda-module.  Qualifying clusters are recorded as the hierarchy grows.
    The final flat partition selects large/high-quality recorded modules first;
    this makes the exported Cytoscape table non-overlapping while preserving
    the HC-PIN hierarchy internally.
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
    for (u, v), merge_ecv in sorted(ecv.items(), key=lambda item: (-item[1], item[0])):
        ru, rv = find(u), find(v)
        if ru == rv:
            continue
        candidate = members[ru] | members[rv]
        if len(members[ru]) < len(members[rv]):
            ru, rv = rv, ru
        parent[rv] = ru
        members[ru] = candidate
        del members[rv]
        if len(candidate) >= min_size and module_condition(candidate, adjacency, mode, threshold):
            internal = [adjacency[x][y] for x in candidate for y in adjacency[x]
                        if x < y and y in candidate]
            candidates.append((set(candidate), merge_ecv,
                               float(np.mean(internal)) if internal else 0.0))

    # Prefer larger modules, then modules born at a stronger ECV and with higher
    # internal similarity. Descendant modules fully covered by a selected parent
    # are omitted from the flat export.
    candidates.sort(key=lambda x: (-len(x[0]), -x[1], -x[2], min(x[0])))
    assignment, module_id, claimed = {u: 0 for u in adjacency}, 0, set()
    for cluster, _, _ in candidates:
        available = cluster - claimed
        if len(available) < min_size:
            continue
        module_id += 1
        for node in available:
            assignment[node] = module_id
        claimed.update(available)
    return assignment, ecv


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


def write_outputs(outdir, spectra, edges, adjacency, matches, assignment, ecv, abundance):
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

    edge_fields = EDGE_FIELDS + ["Edge_Clustering_Value", "Source_Module", "Target_Module", "Same_Module"]
    with open(outdir/"03_edges_with_modules.csv", "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=edge_fields); w.writeheader()
        for u, v, weight, nmatch in edges:
            su, sv = spec_by_id.get(u), spec_by_id.get(v)
            w.writerow({"Source_Scan": u, "Target_Scan": v,
                        "Source_Precursor_mz": su.precursor_mz if su else "",
                        "Target_Precursor_mz": sv.precursor_mz if sv else "",
                        "Precursor_Delta_mz": (sv.precursor_mz-su.precursor_mz) if su and sv else "",
                        "Matched_Peaks": nmatch, "Modified_Cosine_Score": weight,
                        "Edge_Clustering_Value": ecv.get(tuple(sorted((u, v))), ""),
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
        for u, v, score, nmatch in edges:
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
    p.add_argument("--output", default="molecular_network_hcpin_output")
    p.add_argument("--min-cosine", type=float, default=.70)
    p.add_argument("--min-matched-peaks", type=int, default=6)
    p.add_argument("--top-k", type=int, default=5,
                   help="Keep strongest K neighbors per node before HC-PIN; 0 disables")
    p.add_argument("--top-k-mode", choices=["mutual", "union"], default="mutual")
    p.add_argument("--fragment-tolerance", type=float, default=.01)
    p.add_argument("--mode", choices=["weak", "strong"], default="weak")
    p.add_argument("--hcpin-threshold", type=float, default=1.0)
    p.add_argument("--min-module-size", type=int, default=3)
    return p.parse_args()


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
    nodes = {s.scan for s in spectra}
    adjacency, matches = graph_from_edges(nodes, edges)
    print(f"Retained edges: {len(edges)}")
    assignment, ecv = hcpin_partition(adjacency, args.mode, args.hcpin_threshold,
                                      args.min_module_size)
    efficacy = activity_efficacy(args.activity_table, args.activity_group_column,
                                 args.activity_response_column)
    print("Recovery indices:", ", ".join(f"{k}={v:.4f}" for k, v in efficacy.items()))
    abundance = read_abundance_and_label(args.abundance_table, efficacy)
    write_outputs(outdir, spectra, edges, adjacency, matches, assignment, ecv, abundance)
    module_count = len({x for x in assignment.values() if x > 0})
    focused = len({assignment[n] for n in assignment if assignment[n] > 0 and abundance.get(n, {}).get("D_candidate")})
    print(f"HC-PIN modules: {module_count}; modules containing D candidates: {focused}")
    print(f"Output: {outdir.resolve()}")


if __name__ == "__main__":
    main()
