"""Input/output helpers for CIF and POSCAR-like structure files."""
from __future__ import annotations

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]


def ensure_path(path: PathLike) -> Path:
    return Path(path).expanduser().resolve()


def read_atoms(path: PathLike):
    """Read CIF or POSCAR/CONTCAR as ASE Atoms."""
    from ase.io import read

    path = ensure_path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return read(str(path))


def read_structure(path: PathLike):
    """Read CIF or POSCAR/CONTCAR as pymatgen Structure."""
    from pymatgen.core import Structure

    path = ensure_path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return Structure.from_file(str(path))


def atoms_to_structure(atoms):
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_structure(atoms)


def structure_to_atoms(structure):
    from pymatgen.io.ase import AseAtomsAdaptor

    return AseAtomsAdaptor.get_atoms(structure)


def write_atoms(atoms, path: PathLike) -> Path:
    path = ensure_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    atoms.write(str(path))
    return path
