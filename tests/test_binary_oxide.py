from pymatgen.core import Composition

from oxidesse.binary_oxide import balance_binary_oxide_references


def test_llzo_balance():
    refs = balance_binary_oxide_references(Composition("Li7La3Zr2O12"))
    got = {r.formula: r.coefficient for r in refs}
    assert got["Li2O"] == 3.5
    assert got["La2O3"] == 1.5
    assert got["ZrO2"] == 2.0
