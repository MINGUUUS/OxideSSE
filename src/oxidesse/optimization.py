"""Geometry optimization using ASE-compatible calculators."""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .calculators import get_calculator
from .io import read_atoms, write_atoms

LOGGER = logging.getLogger(__name__)


@dataclass
class OptimizationResult:
    atoms: object | None
    converged: bool
    energy: float | None
    energy_per_atom: float | None
    max_force: float | None
    output_path: Path | None
    log_path: Path | None
    message: str = ""


def _max_force(atoms) -> float:
    forces = atoms.get_forces()
    return float(np.sqrt((forces**2).sum(axis=1).max()))


def optimize_atoms(
    atoms,
    calculator="7net-0",
    fmax: float = 0.01,
    max_steps: int = 500,
    cell_relax: bool = True,
    hydrostatic_strain: bool = True,
    fix_symmetry: bool = True,
    optimizer_cls=None,
    log_file: str | Path | None = None,
    device: str = "auto",
    large_force_threshold: float = 1000.0,
) -> OptimizationResult:
    """Optimize an ASE Atoms object.

    Parameters
    ----------
    atoms
        ASE Atoms object.
    calculator
        ASE calculator object or a supported calculator name.
    fmax
        Force convergence criterion in eV/Angstrom.
    max_steps
        Maximum optimizer steps.
    cell_relax
        If True, optimize cell and positions using ``ExpCellFilter``.
    hydrostatic_strain
        Passed to ``ExpCellFilter`` for cell relaxation.
    fix_symmetry
        If True, apply ASE ``FixSymmetry`` before optimization.
    """
    from ase.constraints import FixSymmetry
    from ase.filters import ExpCellFilter
    from ase.optimize import FIRE

    opt_atoms = atoms.copy()
    calc = get_calculator(calculator, device=device)
    opt_atoms.calc = calc

    if fix_symmetry:
        opt_atoms.set_constraint(FixSymmetry(opt_atoms))

    optimizer_cls = optimizer_cls or FIRE
    log_path = Path(log_file).expanduser().resolve() if log_file else None
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if cell_relax:
            target = ExpCellFilter(opt_atoms, hydrostatic_strain=hydrostatic_strain)
            optimizer = optimizer_cls(target, logfile=str(log_path) if log_path else None)
        else:
            optimizer = optimizer_cls(opt_atoms, logfile=str(log_path) if log_path else None)

        converged = bool(optimizer.run(fmax=fmax, steps=max_steps))
        opt_atoms.wrap()
        max_force = _max_force(opt_atoms)
        if max_force > large_force_threshold:
            return OptimizationResult(
                atoms=None,
                converged=False,
                energy=None,
                energy_per_atom=None,
                max_force=max_force,
                output_path=None,
                log_path=log_path,
                message=f"Optimization failed: max force {max_force:.3e} eV/Angstrom exceeds threshold.",
            )
        energy = float(opt_atoms.get_potential_energy())
        return OptimizationResult(
            atoms=opt_atoms,
            converged=converged,
            energy=energy,
            energy_per_atom=energy / len(opt_atoms),
            max_force=max_force,
            output_path=None,
            log_path=log_path,
            message="converged" if converged else "not converged",
        )
    except Exception as exc:
        LOGGER.exception("Optimization failed")
        return OptimizationResult(
            atoms=None,
            converged=False,
            energy=None,
            energy_per_atom=None,
            max_force=None,
            output_path=None,
            log_path=log_path,
            message=str(exc),
        )


def optimize_structure(
    structure: str | Path,
    calculator="7net-0",
    output_dir: str | Path = ".",
    output_filename: str | None = None,
    fmax: float = 0.01,
    max_steps: int = 500,
    cell_relax: bool = True,
    device: str = "auto",
    **kwargs,
) -> OptimizationResult:
    """Read a CIF/POSCAR, optimize it, and write an optimized CIF."""
    structure = Path(structure).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    atoms = read_atoms(structure)
    log_path = output_dir / f"{structure.stem}.opt.log"
    result = optimize_atoms(
        atoms=atoms,
        calculator=calculator,
        fmax=fmax,
        max_steps=max_steps,
        cell_relax=cell_relax,
        log_file=log_path,
        device=device,
        **kwargs,
    )
    if result.atoms is not None:
        output_filename = output_filename or f"opt_{structure.stem}.cif"
        result.output_path = write_atoms(result.atoms, output_dir / output_filename)
    return result


def optimize_structures(
    input_dir: str | Path,
    pattern: str = "*.cif",
    output_dir: str | Path = "optimized",
    calculator="7net-0",
    **kwargs,
) -> list[OptimizationResult]:
    """Batch optimize structures in a directory."""
    input_dir = Path(input_dir).expanduser().resolve()
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files matching {pattern!r} in {input_dir}")
    return [
        optimize_structure(f, calculator=calculator, output_dir=output_dir, **kwargs)
        for f in files
    ]
