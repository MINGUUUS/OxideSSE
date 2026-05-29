"""Binary-oxide-referenced formation energy using MLP recalculated energies."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from .binary_oxide import BinaryOxideReference, balance_binary_oxide_references
from .calculators import calculator_label, get_calculator
from .io import read_structure, structure_to_atoms


@dataclass
class ReferenceEnergy:
    element: str
    formula: str
    coefficient: float
    material_id: str | None
    energy: float
    energy_per_atom: float
    energy_per_formula: float
    nsites: int


@dataclass
class FormationEnergyResult:
    formula: str
    target_energy: float
    target_energy_per_atom: float
    formation_energy: float
    formation_energy_per_atom: float
    references: list[ReferenceEnergy]
    csv_path: Path | None = None


def _get_api_key(api_key: str | None) -> str:
    key = api_key or os.environ.get("MP_API_KEY")
    if not key:
        raise ValueError("Materials Project API key is required. Pass api_key=... or set MP_API_KEY.")
    return key


def _energy_from_structure(structure, calculator) -> float:
    atoms = structure_to_atoms(structure)
    atoms.calc = calculator
    return float(atoms.get_potential_energy())


def _fetch_most_stable_structure(formula: str, api_key: str):
    try:
        from mp_api.client import MPRester
    except ImportError as exc:  # pragma: no cover
        raise ImportError("mp-api is required. Install it with `pip install mp-api`, or reinstall OxideSSE with its default dependencies.") from exc

    with MPRester(api_key) as mpr:
        results = mpr.materials.summary.search(formula=formula, fields=[
            "material_id", "formula_pretty", "energy_above_hull", "structure", "nsites"
        ])
    if not results:
        raise ValueError(f"No Materials Project structure found for formula {formula!r}.")
    best = min(results, key=lambda x: float(x.energy_above_hull or 0.0))
    return best.material_id, best.structure


def compute_binary_oxide_formation_energy(
    structure: str | Path | Any,
    calculator: str | Any = "7net-0",
    api_key: str | None = None,
    output_csv: str | Path | None = None,
    custom_reference_formulas: Mapping[str, str] | None = None,
    device: str = "auto",
    append: bool = True,
) -> FormationEnergyResult:
    """Compute binary-oxide-referenced formation energy for an oxide structure.

    All reference energies, including Li2O, are recalculated with the provided
    ASE-compatible calculator. The input must contain oxygen.
    """
    pmg_structure = read_structure(structure) if isinstance(structure, (str, Path)) else structure
    refs = balance_binary_oxide_references(
        pmg_structure.composition, custom_reference_formulas=custom_reference_formulas
    )
    calc = get_calculator(calculator, device=device)
    api_key = _get_api_key(api_key)

    target_energy = _energy_from_structure(pmg_structure, calc)
    target_energy_per_atom = target_energy / pmg_structure.num_sites

    reference_energies: list[ReferenceEnergy] = []
    ref_sum = 0.0
    from pymatgen.core import Composition

    for ref in refs:
        material_id, ref_structure = _fetch_most_stable_structure(ref.formula, api_key=api_key)
        ref_energy = _energy_from_structure(ref_structure, calc)
        ref_energy_per_atom = ref_energy / ref_structure.num_sites
        atoms_per_ref_formula = Composition(ref.formula).num_atoms
        ref_energy_per_formula = ref_energy_per_atom * atoms_per_ref_formula
        ref_sum += ref.coefficient * ref_energy_per_formula
        reference_energies.append(
            ReferenceEnergy(
                element=ref.element,
                formula=ref.formula,
                coefficient=ref.coefficient,
                material_id=str(material_id),
                energy=ref_energy,
                energy_per_atom=ref_energy_per_atom,
                energy_per_formula=ref_energy_per_formula,
                nsites=ref_structure.num_sites,
            )
        )

    formation_energy = target_energy - ref_sum
    result = FormationEnergyResult(
        formula=pmg_structure.composition.reduced_formula,
        target_energy=target_energy,
        target_energy_per_atom=target_energy_per_atom,
        formation_energy=formation_energy,
        formation_energy_per_atom=formation_energy / pmg_structure.num_sites,
        references=reference_energies,
    )

    if output_csv is not None:
        result.csv_path = _write_formation_csv(result, output_csv, calculator_label(calculator), append=append)
    return result


def _write_formation_csv(
    result: FormationEnergyResult,
    output_csv: str | Path,
    calculator: str,
    append: bool = True,
) -> Path:
    output_csv = Path(output_csv).expanduser().resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    refs = ";".join(
        f"{r.formula}:{r.coefficient:g}:{r.material_id}:E_formula={r.energy_per_formula:.12g}"
        for r in result.references
    )
    row = {
        "formula": result.formula,
        "calculator": calculator,
        "target_energy_eV": result.target_energy,
        "target_energy_per_atom_eV": result.target_energy_per_atom,
        "formation_energy_eV": result.formation_energy,
        "formation_energy_per_atom_eV": result.formation_energy_per_atom,
        "references": refs,
    }
    df = pd.DataFrame([row])
    mode = "a" if append and output_csv.exists() else "w"
    df.to_csv(output_csv, mode=mode, header=(mode == "w"), index=False)
    return output_csv
