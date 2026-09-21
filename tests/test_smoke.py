import csv
import sys
from pathlib import Path

import numpy as np

from ratio_dia.pipeline import (
    activity_efficacy,
    build_relation_aware_graph,
    ecv_consensus_partition,
    hcpin_partition,
    parse_args,
    read_mgf,
    relation_partition,
    Spectrum,
)


ROOT = Path(__file__).resolve().parents[1]


def test_demo_mgf_ids():
    spectra = read_mgf(ROOT / "examples" / "demo.mgf")
    assert len(spectra) == 8
    assert {x.scan for x in spectra} == set(range(1, 9))


def test_recovery_order():
    efficacy = activity_efficacy(ROOT / "examples" / "activity_response.csv")
    assert efficacy["D"] > efficacy["E"] > efficacy["B"] > efficacy["C"] > efficacy["A"]


def _two_dense_groups_with_bridge():
    adjacency = {node: {} for node in range(1, 9)}

    def add(u, v, weight):
        adjacency[u][v] = weight
        adjacency[v][u] = weight

    for group in (range(1, 5), range(5, 9)):
        for u in group:
            for v in group:
                if u < v:
                    add(u, v, 0.95)
    add(4, 5, 0.70)
    return adjacency


def test_quality_cut_splits_a_connected_component():
    adjacency = _two_dense_groups_with_bridge()
    quality, _ = hcpin_partition(adjacency, min_size=3,
                                 selection_mode="quality-cut")
    legacy, _ = hcpin_partition(adjacency, min_size=3,
                                selection_mode="largest-first")
    assert len({value for value in quality.values() if value}) == 2
    assert len({value for value in legacy.values() if value}) == 1
    assert quality[1] == quality[4]
    assert quality[5] == quality[8]
    assert quality[1] != quality[5]


def test_quality_cut_modules_are_disjoint_and_lambda_eligible():
    adjacency = _two_dense_groups_with_bridge()
    assignment, _ = hcpin_partition(adjacency, min_size=3,
                                    selection_mode="quality-cut")
    groups = {}
    for node, module_id in assignment.items():
        if module_id:
            groups.setdefault(module_id, set()).add(node)
    assert set().union(*groups.values()) == set(adjacency)
    assert sum(map(len, groups.values())) == len(adjacency)
    from ratio_dia.pipeline import module_condition
    assert all(module_condition(group, adjacency, "weak", 1.0)
               for group in groups.values())


def test_primary_cli_disables_top_k_by_default(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "ratio-dia",
        "--mgf", "study.mgf",
        "--abundance-table", "abundance.csv",
        "--activity-table", "activity.csv",
    ])
    assert parse_args().top_k == 0
    assert parse_args().partition_method == "relation-aware"
    assert parse_args().input_scope == "edge-connected"
    assert parse_args().random_seed == 20260914


def test_ecv_consensus_is_deterministic_and_retains_edge_evidence():
    adjacency = _two_dense_groups_with_bridge()
    first, ecv, evidence = ecv_consensus_partition(
        adjacency, consensus_runs=5, final_seed=123
    )
    second, _, _ = ecv_consensus_partition(
        adjacency, consensus_runs=5, final_seed=123
    )
    assert first == second
    assert set(evidence) == set(ecv)
    assert all(0 <= row["Consensus_Coassignment"] <= 1
               for row in evidence.values())


def test_relation_graph_adds_nonreciprocal_motif_only_relationship():
    def spectrum(scan, precursor, peaks):
        mz = np.asarray(peaks, dtype=float)
        intensity = np.ones(len(mz), dtype=float)
        intensity /= np.linalg.norm(intensity)
        return Spectrum(scan, precursor, 1.0, mz, intensity)

    spectra = [
        spectrum(1, 300, [60, 100, 150]),
        spectrum(2, 301, [60, 100, 151]),
        spectrum(3, 302, [70, 120, 160]),
        spectrum(4, 303, [80, 130, 170]),
    ]
    original = [(1, 3, .8, 6), (2, 4, .8, 6), (3, 4, .8, 6)]
    graph, _, _, _ = build_relation_aware_graph(
        {1, 2, 3, 4}, spectra, original, motif_k=1,
        motif_minimum=.05, pair_motif_weight=0,
    )
    assert graph.has_edge(1, 2)
    assert graph[1][2]["Edge_Origin"] == "motif_neighborhood"
    assert graph[1][2]["Modified_Cosine_Score"] == 0


def test_relation_partition_is_deterministic():
    import networkx as nx
    graph = nx.Graph()
    graph.add_weighted_edges_from([
        (1, 2, 2.0), (1, 3, 2.0), (2, 3, 2.0),
        (4, 5, 2.0), (4, 6, 2.0), (5, 6, 2.0),
        (3, 4, 0.05),
    ])
    _, first = relation_partition(graph, resolution=1.0, seed=77)
    _, second = relation_partition(graph, resolution=1.0, seed=77)
    assert first == second
