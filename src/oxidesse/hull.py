"""MLP-recalculated phase diagram and energy-above-hull calculations."""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .calculators import calculator_label, get_calculator
from .formation_energy import _get_api_key
from .io import read_structure, structure_to_atoms


@dataclass
class HullEnergyResult:
    formula: str
    chemical_system: str
    thermo_type: str | None
    target_energy: float
    target_energy_per_atom: float
    energy_above_hull: float
    n_entries: int
    failed_entries: list[str]
    csv_path: Path | None = None
    cache_path: Path | None = None


def _structure_energy(structure, calculator) -> float:
    atoms = structure_to_atoms(structure)
    atoms.calc = calculator
    return float(atoms.get_potential_energy())


def _entry_key(entry) -> str:
    return str(getattr(entry, "entry_id", None) or entry.data.get("material_id", None) or entry.composition.reduced_formula)


def _load_cache(path: Path) -> dict:
    if path.exists():
        return json.loads(path.read_text())
    return {}


def _save_cache(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True))


def _normalize_thermo_type(thermo_type: str | None) -> str | None:
    """Validate and normalize Materials Project thermo_type names."""
    if thermo_type is None:
        return None
    allowed = {"GGA_GGA+U", "R2SCAN", "GGA_GGA+U_R2SCAN"}
    value = str(thermo_type).upper()
    aliases = {
        "GGA": "GGA_GGA+U",
        "GGAU": "GGA_GGA+U",
        "GGA+U": "GGA_GGA+U",
        "GGA_GGAU": "GGA_GGA+U",
        "GGA_GGA_U": "GGA_GGA+U",
        "GGA_GGA+U": "GGA_GGA+U",
        "R2SCAN": "R2SCAN",
        "GGAU_R2SCAN": "GGA_GGA+U_R2SCAN",
        "GGA_GGAU_R2SCAN": "GGA_GGA+U_R2SCAN",
        "GGA_GGA+U_R2SCAN": "GGA_GGA+U_R2SCAN",
        "GGA_R2SCAN": "GGA_GGA+U_R2SCAN",
    }
    normalized = aliases.get(value, thermo_type)
    if normalized not in allowed:
        raise ValueError(
            f"Unsupported thermo_type={thermo_type!r}. "
            f"Use one of {sorted(allowed)} or None."
        )
    return normalized


def _thermo_cache_label(thermo_type: str | None) -> str:
    """Return a short filesystem-safe label for a normalized thermo_type."""
    if thermo_type is None:
        return "all"
    labels = {
        "GGA_GGA+U": "GGAU",
        "R2SCAN": "R2SCAN",
        "GGA_GGA+U_R2SCAN": "GGAU_R2SCAN",
    }
    return labels.get(thermo_type, str(thermo_type).replace("+", "U"))


def _fetch_mp_entries(chemsys: list[str], api_key: str, thermo_type: str | None = "GGA_GGA+U"):
    """Fetch MP entries for a chemical system with an explicit thermo_type filter.

    Use mp_api.client.MPRester directly. Older/legacy pymatgen MPRester
    wrappers may ignore or fail to support additional_criteria, which can make
    thermo_type changes appear to have no effect.
    """
    try:
        from mp_api.client import MPRester
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Materials Project support is required. Install mp-api, or install "
            "OxideSSE with its default dependencies."
        ) from exc

    thermo_type = _normalize_thermo_type(thermo_type)
    criteria = {"thermo_types": [thermo_type]} if thermo_type else None
    with MPRester(api_key) as mpr:
        # inc_structure=True is important because OxideSSE replaces every MP
        # entry energy with an MLP energy evaluated on the MP structure.
        kwargs = {"inc_structure": True}
        if criteria is not None:
            kwargs["additional_criteria"] = criteria
        return mpr.get_entries_in_chemsys(chemsys, **kwargs)


def _ensure_entry_structure(entry, api_key: str):
    if getattr(entry, "structure", None) is not None:
        return entry.structure
    try:
        from mp_api.client import MPRester
    except ImportError as exc:  # pragma: no cover
        raise ImportError("mp-api is required to fetch structures for entries without structures.") from exc

    mid = _entry_key(entry)
    with MPRester(api_key) as mpr:
        doc = mpr.materials.summary.search(material_ids=[mid], fields=["structure"])
    if not doc:
        raise ValueError(f"Could not fetch structure for entry {mid}.")
    return doc[0].structure


def compute_energy_above_hull(
    structure: str | Path | Any,
    calculator: str | Any = "7net-0",
    api_key: str | None = None,
    thermo_type: str | None = "GGA_GGA+U",
    output_csv: str | Path | None = None,
    cache_dir: str | Path | None = ".oxidesse_cache",
    use_cache: bool = True,
    device: str = "auto",
    append: bool = True,
) -> HullEnergyResult:
    """Compute energy above hull using MP structures but MLP-recalculated energies.

    MP entries are used to define the chemical system and phase diagram entries.
    Their energies are replaced by energies evaluated with the provided MLP/ASE
    calculator before building the phase diagram.
    """
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from pymatgen.entries.computed_entries import ComputedEntry

    pmg_structure = read_structure(structure) if isinstance(structure, (str, Path)) else structure
    calc = get_calculator(calculator, device=device)
    calc_label = calculator_label(calculator)
    api_key = _get_api_key(api_key)
    thermo_type = _normalize_thermo_type(thermo_type)

    formula = pmg_structure.composition.reduced_formula
    chemsys_list = sorted([el.symbol for el in pmg_structure.composition.elements])
    chemical_system = "-".join(chemsys_list)

    target_energy = _structure_energy(pmg_structure, calc)
    target_entry = ComputedEntry(pmg_structure.composition, target_energy, entry_id="target")

    cache_path = None
    cache = {}
    if cache_dir is not None:
        safe_thermo = _thermo_cache_label(thermo_type)
        cache_path = Path(cache_dir).expanduser().resolve() / f"{chemical_system}_{calc_label}_{safe_thermo}.json"
        cache = _load_cache(cache_path) if use_cache else {}

    mp_entries = _fetch_mp_entries(chemsys_list, api_key=api_key, thermo_type=thermo_type)
    if not mp_entries:
        raise ValueError(
            f"No Materials Project entries were fetched for chemical system {chemical_system} "
            f"with thermo_type={thermo_type!r}. Try different thermo_type, check your API key and network connection."
        )

    updated_entries = []
    failed_entries: list[str] = []

    for entry in mp_entries:
        key = _entry_key(entry)
        cache_key = f"{calc_label}|{key}"
        try:
            if use_cache and cache_key in cache:
                energy = float(cache[cache_key]["energy_eV"])
            else:
                entry_structure = _ensure_entry_structure(entry, api_key=api_key)
                energy = _structure_energy(entry_structure, calc)
                if cache_path is not None:
                    cache[cache_key] = {
                        "entry_id": key,
                        "formula": entry.composition.reduced_formula,
                        "energy_eV": energy,
                        "energy_per_atom_eV": energy / entry.composition.num_atoms,
                    }
            updated_entries.append(ComputedEntry(entry.composition, energy, entry_id=key))
        except Exception as exc:
            failed_entries.append(f"{key}: {exc}")

    if cache_path is not None:
        _save_cache(cache_path, cache)

    if not updated_entries:
        details = "\n".join(failed_entries[:10])
        raise ValueError(
            "Materials Project entries were fetched, but none could be recalculated "
            "with the MLP calculator. First failed entries:\n" + details
        )

    phase_diagram = PhaseDiagram(updated_entries + [target_entry])
    energy_above_hull = float(phase_diagram.get_e_above_hull(target_entry))

    result = HullEnergyResult(
        formula=formula,
        chemical_system=chemical_system,
        thermo_type=thermo_type,
        target_energy=target_energy,
        target_energy_per_atom=target_energy / pmg_structure.num_sites,
        energy_above_hull=energy_above_hull,
        n_entries=len(updated_entries),
        failed_entries=failed_entries,
        cache_path=cache_path,
    )
    if output_csv is not None:
        result.csv_path = _write_hull_csv(result, output_csv, calc_label, append=append)
    return result


def _write_hull_csv(result: HullEnergyResult, output_csv: str | Path, calculator: str, append: bool = True) -> Path:
    output_csv = Path(output_csv).expanduser().resolve()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "formula": result.formula,
        "chemical_system": result.chemical_system,
        "thermo_type": result.thermo_type or "all",
        "calculator": calculator,
        "target_energy_eV": result.target_energy,
        "target_energy_per_atom_eV": result.target_energy_per_atom,
        "energy_above_hull_eV_per_atom": result.energy_above_hull,
        "n_entries": result.n_entries,
        "n_failed_entries": len(result.failed_entries),
        "failed_entries": ";".join(result.failed_entries),
        "cache_path": str(result.cache_path) if result.cache_path else "",
    }
    df = pd.DataFrame([row])
    mode = "a" if append and output_csv.exists() else "w"
    df.to_csv(output_csv, mode=mode, header=(mode == "w"), index=False)
    return output_csv
