"""Small CLI wrapper for OxideSSE."""
from __future__ import annotations

from .arrhenius import plot_arrhenius
from .diffusion import compute_diffusivity_from_lammps
from .formation_energy import compute_binary_oxide_formation_energy
from .hull import compute_energy_above_hull_mlp
from .optimization import optimize_structure, optimize_structures


def main():
    try:
        from fire import Fire
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "The OxideSSE command-line interface requires the 'fire' dependency. "
            "Install it with `pip install fire`, or reinstall OxideSSE with its default dependencies."
        ) from exc

    Fire({
        "optimize-structure": optimize_structure,
        "optimize-structures": optimize_structures,
        "formation-energy": compute_binary_oxide_formation_energy,
        "energy-above-hull": compute_energy_above_hull_mlp,
        "diffusivity": compute_diffusivity_from_lammps,
        "arrhenius": plot_arrhenius,
    })


if __name__ == "__main__":
    main()
