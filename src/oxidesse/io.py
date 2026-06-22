"""Input/output helpers for CIF, POSCAR/CONTCAR, and LAMMPS data files."""
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


def read_lammps_data_structure(path: PathLike, atom_style: str = "atomic"):
    """Read a LAMMPS data file as a pymatgen Structure."""
    from pymatgen.io.lammps.data import LammpsData

    path = ensure_path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return LammpsData.from_file(str(path), atom_style=atom_style).structure


def read_structure(path: PathLike):
    """Read CIF, POSCAR/CONTCAR, or LAMMPS data as pymatgen Structure.

    CIF and POSCAR/CONTCAR-like files are first parsed with pymatgen's
    Structure.from_file.  Files with common LAMMPS data extensions, or files
    that fail the standard structure parser, are parsed as LAMMPS data files
    using atom_style="atomic".
    """
    from pymatgen.core import Structure

    path = ensure_path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    lower_name = path.name.lower()
    lammps_like = lower_name.endswith((".data", ".lmp", ".lammps"))
    if lammps_like:
        return read_lammps_data_structure(path)

    try:
        return Structure.from_file(str(path))
    except Exception as structure_exc:
        try:
            return read_lammps_data_structure(path)
        except Exception as lammps_exc:
            raise ValueError(
                f"Could not read {path} as CIF/POSCAR/CONTCAR or LAMMPS data."
            ) from lammps_exc


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
