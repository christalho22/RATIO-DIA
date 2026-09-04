"""Evaluate chemical coherence of topology-focused modules using all candidates."""
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean, median

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem
    from rdkit.Chem.Scaffolds import MurckoScaffold
except ImportError as exc:  # pragma: no cover
    raise SystemExit("Install the optional dependency with: pip install -e '.[structure]'") from exc


def allowed_formula(formula: str) -> bool:
    elements = set(re.findall(r"[A-Z][a-z]?", formula or ""))
    return bool(elements) and elements <= {"C", "H", "O", "N"}


def load_candidates(path: Path):
    result = defaultdict(list)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if not allowed_formula(row.get("Formula", "")):
                continue
            mol = Chem.MolFromSmiles(row.get("SMILES", ""))
            if mol is None:
                continue
            rec = dict(row)
            rec["mol"] = mol
            rec["fp"] = AllChem.GetMorganGenerator(radius=2, fpSize=2048).GetFingerprint(mol)
            rec["scaffold"] = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or Chem.MolToSmiles(mol)
            result[str(int(float(row["Alignment_ID"])))].append(rec)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", required=True, help="RATIO-DIA node CSV containing Alignment_ID and Module_ID")
    parser.add_argument("--candidates", required=True, help="Long candidate CSV; see docs/input-formats.md")
    parser.add_argument("--output", default="structure_consistency.csv")
    parser.add_argument("--min-module-size", type=int, default=11)
    args = parser.parse_args()

    with open(args.nodes, encoding="utf-8-sig", newline="") as handle:
        nodes = list(csv.DictReader(handle))
    modules = defaultdict(list)
    for row in nodes:
        mid = int(float(row["Module_ID"]))
        if mid > 0:
            modules[mid].append(str(int(float(row["Alignment_ID"]))))
    candidates = load_candidates(Path(args.candidates))
    rows = []
    for mid, ids in sorted(modules.items()):
        if len(ids) < args.min_module_size:
            continue
        annotated = [node for node in ids if candidates.get(node)]
        similarities, same_scaffold, same_ontology = [], 0, 0
        for left, right in combinations(annotated, 2):
            ca, cb = candidates[left], candidates[right]
            best = max(DataStructs.TanimotoSimilarity(a["fp"], b["fp"]) for a in ca for b in cb)
            similarities.append(best)
            same_scaffold += any(a["scaffold"] == b["scaffold"] for a in ca for b in cb)
            same_ontology += any(a.get("Ontology") and a.get("Ontology") == b.get("Ontology")
                                 for a in ca for b in cb)
        pairs = len(similarities)
        ontology_counts = Counter(c.get("Ontology", "") for node in annotated for c in candidates[node]
                                  if c.get("Ontology"))
        rows.append({
            "Module_ID": mid, "Node_Count": len(ids), "Annotated_CHON_Nodes": len(annotated),
            "Annotation_Coverage": len(annotated) / len(ids) if ids else 0,
            "Evaluated_Node_Pairs": pairs,
            "Mean_Max_Tanimoto": mean(similarities) if similarities else "",
            "Median_Max_Tanimoto": median(similarities) if similarities else "",
            "Pairs_MaxTanimoto_ge_0.5": sum(x >= 0.5 for x in similarities) / pairs if pairs else 0,
            "Pairs_MaxTanimoto_ge_0.7": sum(x >= 0.7 for x in similarities) / pairs if pairs else 0,
            "Pairs_AnySame_Murcko": same_scaffold / pairs if pairs else 0,
            "Pairs_AnySame_Ontology": same_ontology / pairs if pairs else 0,
            "Dominant_Candidate_Ontology": ontology_counts.most_common(1)[0][0] if ontology_counts else "",
        })
    if not rows:
        raise SystemExit("No eligible modules were found")
    with open(args.output, "w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(Path(args.output).resolve())


if __name__ == "__main__":
    main()
