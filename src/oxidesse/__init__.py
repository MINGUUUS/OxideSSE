"""OxideSSE: tools for oxide solid-state electrolyte screening."""

__version__ = "0.1.18"


from .arrhenius import ArrheniusResult, fit_arrhenius, plot_arrhenius
from .diffusion import DiffusivityResult, compute_diffusivity_from_lammps
from .formation_energy import FormationEnergyResult, compute_binary_oxide_formation_energy
from .hull import HullEnergyResult, compute_energy_above_hull
from .optimization import OptimizationResult, optimize_structure, optimize_structures

__all__ = [
    "ArrheniusResult",
    "DiffusivityResult",
    "FormationEnergyResult",
    "HullEnergyResult",
    "OptimizationResult",
    "compute_binary_oxide_formation_energy",
    "compute_diffusivity_from_lammps",
    "compute_energy_above_hull",
    "fit_arrhenius",
    "optimize_structure",
    "optimize_structures",
    "plot_arrhenius",
    "__version__",
]
