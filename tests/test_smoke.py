import csv
from pathlib import Path

from ratio_dia.pipeline import activity_efficacy, read_mgf


ROOT = Path(__file__).resolve().parents[1]


def test_demo_mgf_ids():
    spectra = read_mgf(ROOT / "examples" / "demo.mgf")
    assert len(spectra) == 8
    assert {x.scan for x in spectra} == set(range(1, 9))


def test_recovery_order():
    efficacy = activity_efficacy(ROOT / "examples" / "activity_response.csv")
    assert efficacy["D"] > efficacy["E"] > efficacy["B"] > efficacy["C"] > efficacy["A"]
