"""Binary oxide reference selection and composition balancing."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

B3_ELEMENTS = [
    "La", "B", "N", "Al", "P", "Ga", "As", "Br", "In", "Sb", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Tl",
    "Bi", "Ac", "Sc", "Cr", "Fe", "Co", "Y", "Rh", "Au",
]

C4_ELEMENTS = [
    "Si", "S", "Ge", "Se", "Sn", "Te", "Ce", "Pb", "Th", "Zr", "Ti", "Mn",
    "Mo", "Tc", "Ru", "Pd", "Hf", "W", "Re", "Os", "Ir", "Pt",
]


@dataclass(frozen=True)
class BinaryOxideReference:
    element: str
    formula: str
    coefficient: float


def reference_formula_for_element(element: str) -> str:
    """Return the target binary oxide formula for a non-oxygen element."""
    if element == "O":
        raise ValueError("Oxygen is not represented by a separate binary oxide reference.")
    if element == "Li":
        return "Li2O"
    if element in B3_ELEMENTS:
        return f"{element}2O3"
    if element in C4_ELEMENTS:
        return f"{element}O2"
    raise ValueError(
        f"No binary oxide rule is defined for element {element!r}. Add it to "
        "B3_ELEMENTS/C4_ELEMENTS or pass custom reference formulas."
    )


def balance_binary_oxide_references(
    composition,
    custom_reference_formulas: Mapping[str, str] | None = None,
    oxygen_tolerance: float = 1e-6,
) -> list[BinaryOxideReference]:
    """Balance an oxide composition against binary oxide references.

    The coefficient for each reference is set by matching the amount of the
    corresponding non-oxygen element in the target composition. The oxygen count
    is then checked. If the binary oxide references do not reproduce the target
    oxygen amount, a ValueError is raised.
    """
    from pymatgen.core import Composition

    comp = Composition(composition)
    amounts = {el.symbol: float(n) for el, n in comp.items()}
    if "O" not in amounts or amounts["O"] <= 0:
        raise ValueError("Input structure must be an oxide and contain oxygen.")

    custom_reference_formulas = dict(custom_reference_formulas or {})
    refs: list[BinaryOxideReference] = []
    oxygen_from_refs = 0.0

    for element, amount in sorted(amounts.items()):
        if element == "O":
            continue
        formula = custom_reference_formulas.get(element) or reference_formula_for_element(element)
        ref_comp = Composition(formula)
        element_amount = float(ref_comp.get_el_amt_dict().get(element, 0.0))
        oxygen_amount = float(ref_comp.get_el_amt_dict().get("O", 0.0))
        if element_amount <= 0 or oxygen_amount <= 0:
            raise ValueError(f"Reference formula {formula!r} is not a valid {element}-oxide.")
        coeff = amount / element_amount
        oxygen_from_refs += coeff * oxygen_amount
        refs.append(BinaryOxideReference(element=element, formula=formula, coefficient=coeff))

    target_oxygen = amounts["O"]
    if abs(oxygen_from_refs - target_oxygen) > oxygen_tolerance:
        raise ValueError(
            "Binary oxide references do not oxygen-balance the target composition: "
            f"target O={target_oxygen:.8g}, references O={oxygen_from_refs:.8g}. "
            "Provide custom_reference_formulas or check oxidation-state assumptions."
        )
    return refs
