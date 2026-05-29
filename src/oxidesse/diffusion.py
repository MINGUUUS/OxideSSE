"""MSD and diffusivity analysis from LAMMPS dump trajectories."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass
class SpeciesDiffusivity:
    species: str
    diffusivity: float
    diffusivity_std: float | None = None
    diffusivity_relative_std: float | None = None
    diffusivity_components: list[float] | None = None
    n_jump: float | None = None
    n_jump_component: list[int] | list[float] | None = None
    fit_start_time_ps: float | None = None
    fit_end_time_ps: float | None = None
    msd_start: float | None = None
    msd_middle: float | None = None
    msd_last: float | None = None


@dataclass
class DiffusivityResult:
    temperature: float
    timestep_fs: float
    step_skip: int
    primary_species: str
    species_results: list[SpeciesDiffusivity]
    msd_csv_path: Path | None = None
    summary_csv_path: Path | None = None
    plot_path: Path | None = None

    @property
    def diffusivity(self) -> float | None:
        for item in self.species_results:
            if item.species == self.primary_species:
                return item.diffusivity
        return None

    @property
    def diffusivity_std(self) -> float | None:
        for item in self.species_results:
            if item.species == self.primary_species:
                return item.diffusivity_std
        return None


def _get_diffusion_classes():
    """Return OxideSSE's vendored AIMD-style diffusion-analysis classes."""
    try:
        from ._aimd import DiffusivityAnalyzer, ErrorAnalysisFromDiffusivityAnalyzer
    except Exception as exc:
        raise ImportError(
            "OxideSSE uses its internal `oxidesse._aimd` diffusion backend; "
            "the external `aimd` package is not required. Please reinstall OxideSSE "
            "and restart the Python/Jupyter kernel if this import fails."
        ) from exc

    return DiffusivityAnalyzer, ErrorAnalysisFromDiffusivityAnalyzer

def _element_from_mass(mass: float, tolerance: float = 0.1) -> str:
    from pymatgen.core.periodic_table import Element

    for el in Element:
        try:
            if abs(float(el.atomic_mass) - float(mass)) <= tolerance:
                return el.symbol
        except TypeError:
            continue
    raise ValueError(f"Could not infer element from atomic mass {mass}.")


def lammps_dump_to_structures(
    simulation_dir: str | Path,
    dump_file: str | Path,
    data_file: str | Path,
    save_poscars: bool = False,
    output_poscar_dir: str | Path | None = None,
) -> tuple[list[Any], list[int], Path | None]:
    """Convert a single LAMMPS dump trajectory to pymatgen Structures.

    Unlike the research notebook, this function does not infer filenames and does
    not apply any supercell transformation.
    """
    from pymatgen.core import Structure
    from pymatgen.io.lammps.data import LammpsData
    from pymatgen.io.lammps.outputs import parse_lammps_dumps

    simulation_dir = Path(simulation_dir).expanduser().resolve()
    dump_path = simulation_dir / dump_file
    data_path = simulation_dir / data_file
    if not dump_path.exists():
        raise FileNotFoundError(dump_path)
    if not data_path.exists():
        raise FileNotFoundError(data_path)

    lammps_data = LammpsData.from_file(str(data_path), atom_style="atomic")
    base_structure = lammps_data.structure
    masses = np.array(lammps_data.masses).flatten().tolist()
    species_by_type = [_element_from_mass(m) for m in masses]

    poscar_dir = None
    if save_poscars:
        poscar_dir = Path(output_poscar_dir or (simulation_dir / "traj_to_POSCARs")).expanduser().resolve()
        poscar_dir.mkdir(parents=True, exist_ok=True)

    structures = []
    timesteps = []
    for dump in parse_lammps_dumps(str(dump_path)):
        frame = dump.data.sort_values("type").copy()
        species = [species_by_type[int(t) - 1] for t in frame["type"].to_numpy()]
        coords = frame[["x", "y", "z"]].to_numpy(dtype=float)
        if len(coords) != len(base_structure):
            raise ValueError("LAMMPS data file and dump frame have different numbers of atoms.")
        struct = Structure(base_structure.lattice, species, coords, coords_are_cartesian=True, to_unit_cell=False)
        structures.append(struct)
        timesteps.append(int(dump.timestep))
        if save_poscars and poscar_dir is not None:
            struct.to(fmt="poscar", filename=str(poscar_dir / f"POSCAR-{int(dump.timestep):07d}"))
    return structures, timesteps, poscar_dir


# AIMD fitting controls used by default for Li-ion diffusion.
# For Li-ion diffusion, early-time ballistic motion should be excluded and a
# sufficiently large MSD should be reached to identify meaningful hopping events.
# See: He, X., Zhu, Y., Epstein, A. et al. Statistical variances of diffusional
# properties from ab initio molecular dynamics simulations. npj Computational
# Materials 4, 18 (2018).
DEFAULT_LI_SPEC_DICT = {
    "lower_bound": 4.5,
    "upper_bound": 0.99,
    "minimum_msd_diff": 4.5,
    "total_sim_time_limit": 5,
}

# Framework/cation/anion species excluding Li should ideally remain immobile.
# They should have very small MSD, which can be used as an indicator for melting
# or structural instability. A separate spec_dict is used to compute and track
# diffusivities of these nearly immobile ions without requiring Li-scale MSD.
DEFAULT_NON_LI_SPEC_DICT = {
    "lower_bound": 0.1,
    "upper_bound": 0.99,
    "minimum_msd_diff": 0.01,
    "total_sim_time_limit": 5,
}


def _default_spec_dict_for_species(species: str) -> dict[str, float]:
    return DEFAULT_LI_SPEC_DICT.copy() if species == "Li" else DEFAULT_NON_LI_SPEC_DICT.copy()


def _complete_spec_dict(species: str, spec_dict: dict[str, float] | None) -> dict[str, float]:
    """Merge user-provided fitting controls with species-specific defaults."""
    base = _default_spec_dict_for_species(species)
    if spec_dict:
        base.update(spec_dict)
    return base


def _analyze_species(
    structures,
    species: str,
    temperature: float,
    timestep_fs: float,
    step_skip: int,
    spec_dict: dict[str, float] | None = None,
):
    DiffusivityAnalyzer, _ = _get_diffusion_classes()
    return DiffusivityAnalyzer.from_structures(
        structures=structures,
        specie=species,
        temperature=temperature,
        time_step=timestep_fs,
        step_skip=step_skip,
        time_intervals_number=len(structures),
        spec_dict=_complete_spec_dict(species, spec_dict),
    )


def _get_aimd_error_summary(difs, oxidized_specie: str) -> dict[str, Any]:
    """Compute AIMD error-analysis summary for one species.

    The original notebooks used
    ``ErrorAnalysisFromDiffusivityAnalyzer(difs).get_summary_dict(...)``
    to obtain jump statistics and diffusivity uncertainties. These values are
    part of the public OxideSSE output because they define the error range used
    in later Arrhenius analysis.
    """
    _, ErrorAnalysisFromDiffusivityAnalyzer = _get_diffusion_classes()
    analyzer = ErrorAnalysisFromDiffusivityAnalyzer(difs)
    return analyzer.get_summary_dict(oxidized_specie=oxidized_specie)


def _as_python_list(value):
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return value


def _as_float_or_none(value):
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def compute_diffusivity_from_lammps(
    simulation_dir: str | Path,
    dump_file: str | Path,
    data_file: str | Path,
    timestep_fs: float,
    step_skip: int,
    temperature: float,
    species: str = "Li",
    output_dir: str | Path = ".",
    save_poscars: bool = False,
    save_plot: bool = True,
    output_csv: str | Path | None = "diffusivity_summary.csv",
    msd_csv: str | Path | None = "msd_by_species.csv",
    plot_other_species: bool = True,
    oxidized_species: dict[str, str] | None = None,
    spec_dict: dict[str, float] | None = None,
    species_spec_dicts: dict[str, dict[str, float]] | None = None,
) -> DiffusivityResult:
    """Compute MSD and diffusivity from one LAMMPS simulation folder.

    ``timestep_fs`` must be in femtoseconds. MSD plots are always reported with
    time in picoseconds.

    Parameters
    ----------
    oxidized_species
        Optional mapping used for conductivity/error analysis, e.g. ``{"Li": "Li+"}``.
        If omitted, only the primary species Li is treated as ``Li+``. Non-Li
        species are analyzed with neutral element symbols and no conductivity conversion.
    spec_dict
        Optional diffusion fitting controls applied to every analyzed species. If omitted,
        OxideSSE uses element-specific defaults: larger MSD thresholds for mobile Li
        and smaller MSD thresholds for nearly immobile non-Li species. The
        default total_sim_time_limit is 5 ps.
    species_spec_dicts
        Optional per-species diffusion fitting controls, e.g. ``{"Li": {...}, "O": {...}}``.
        Values here override both ``spec_dict`` and the built-in defaults.
    """
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    structures, timesteps, poscar_dir = lammps_dump_to_structures(
        simulation_dir=simulation_dir,
        dump_file=dump_file,
        data_file=data_file,
        save_poscars=save_poscars,
        output_poscar_dir=output_dir / "traj_to_POSCARs",
    )
    if len(structures) < 2:
        raise ValueError("At least two trajectory frames are required for diffusivity analysis.")

    elements = sorted({site.specie.symbol for site in structures[0]})
    species_to_analyze = [species] + [el for el in elements if el != species]

    species_results: list[SpeciesDiffusivity] = []
    msd_rows = []
    fig, ax = plt.subplots(figsize=(5, 4)) if save_plot else (None, None)

    for sp in species_to_analyze:
        try:
            current_spec_dict = None
            if species_spec_dicts and sp in species_spec_dicts:
                current_spec_dict = species_spec_dicts[sp]
            elif spec_dict is not None:
                current_spec_dict = spec_dict
            difs = _analyze_species(structures, sp, temperature, timestep_fs, step_skip, spec_dict=current_spec_dict)
            time_ps = np.asarray(difs.dt, dtype=float) / 1000.0
            msd = np.asarray(difs.msd, dtype=float)
            for t, m in zip(time_ps, msd):
                msd_rows.append({"time_ps": t, "species": sp, "msd_A2": m})
            mid = len(msd) // 2
            oxidized_species = oxidized_species or {}
            # The diffusion analyzer itself should always receive the neutral element
            # symbol (e.g., "Li", "Ce", "La", "O") because the trajectory structures
            # are usually not oxidation-state decorated. Oxidized species strings are
            # only needed when converting Li diffusivity to ionic conductivity via the
            # Nernst-Einstein relation. Passing neutral non-Li symbols such as "Ce" or
            # "O" to Specie.from_str() would raise "Invalid species string", so non-Li
            # species default to None unless the user explicitly provides an oxidation
            # state such as {"Ce": "Ce4+"}.
            oxidized_specie = oxidized_species.get(sp)
            if oxidized_specie is None and sp == species and sp == "Li":
                oxidized_specie = "Li+"
            error_summary = _get_aimd_error_summary(difs, oxidized_specie=oxidized_specie)

            diffusivity_std = error_summary.get("diffusivity_standard_deviation")
            diffusivity_relative_std = error_summary.get("diffusivity_relative_standard_deviation")
            n_jump = error_summary.get("n_jump")
            n_jump_component = error_summary.get("n_jump_component")
            components = getattr(difs, "diffusivity_components", None)
            lower_idx = getattr(difs, "lower_bound_index", None)
            upper_idx = getattr(difs, "upper_bound_index", None)
            fit_start = float(time_ps[lower_idx]) if lower_idx is not None and lower_idx < len(time_ps) else None
            fit_end = float(time_ps[upper_idx]) if upper_idx is not None and upper_idx < len(time_ps) else None
            species_results.append(
                SpeciesDiffusivity(
                    species=sp,
                    diffusivity=float(difs.diffusivity),
                    diffusivity_std=_as_float_or_none(diffusivity_std),
                    diffusivity_relative_std=_as_float_or_none(diffusivity_relative_std),
                    diffusivity_components=None if components is None else list(np.asarray(components, dtype=float)),
                    n_jump=_as_float_or_none(n_jump),
                    n_jump_component=_as_python_list(n_jump_component),
                    fit_start_time_ps=fit_start,
                    fit_end_time_ps=fit_end,
                    msd_start=float(msd[0]),
                    msd_middle=float(msd[mid]),
                    msd_last=float(msd[-1]),
                )
            )
            if save_plot and (sp == species or plot_other_species):
                alpha = 1.0 if sp == species else 0.35
                size = 16 if sp == species else 8
                ax.scatter(time_ps, msd, s=size, alpha=alpha, label=sp)

                # Draw the linear fitting region for Li. The AIMD analyzer selects
                # lower/upper indices based on the Li MSD threshold; visualizing this
                # line makes the diffusivity fit region explicit in the saved plot.
                if sp == "Li" and lower_idx is not None and upper_idx is not None and float(difs.diffusivity) > 0:
                    x_fit = np.asarray(time_ps[lower_idx:upper_idx], dtype=float)
                    y_fit = np.asarray(msd[lower_idx:upper_idx], dtype=float)
                    if len(x_fit) >= 2:
                        coeff = np.polyfit(x_fit, y_fit, 1)
                        ax.plot(x_fit, coeff[0] * x_fit + coeff[1], linewidth=2, linestyle="-", alpha=0.8, label="Li fit")
        except Exception as exc:
            species_results.append(SpeciesDiffusivity(species=sp, diffusivity=float("nan")))
            print(f"[{sp}] diffusivity analysis skipped: {exc}")

    plot_path = None
    if save_plot:
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(r"MSD ($\AA^2$)")
        primary_diff = None
        for item in species_results:
            if item.species == species:
                primary_diff = item.diffusivity
                break
        if primary_diff is not None and np.isfinite(primary_diff):
            ax.set_title(f"MSD at {temperature:g} K\n{species} D = {primary_diff:.3e} cm$^2$/s")
        else:
            ax.set_title(f"MSD at {temperature:g} K")
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(frameon=False)
        fig.tight_layout()
        plot_path = output_dir / f"msd_{species}_{temperature:g}K.png"
        fig.savefig(plot_path, dpi=300)
        plt.close(fig)

    msd_csv_path = None
    if msd_csv is not None:
        msd_csv_path = output_dir / msd_csv
        pd.DataFrame(msd_rows).to_csv(msd_csv_path, index=False)

    summary_csv_path = None
    if output_csv is not None:
        summary_csv_path = output_dir / output_csv
        rows = []
        for r in species_results:
            row = {
                "simulation_dir": str(Path(simulation_dir).expanduser().resolve()),
                "dump_file": str(dump_file),
                "data_file": str(data_file),
                "temperature_K": temperature,
                "species": r.species,
                "timestep_fs": timestep_fs,
                "step_skip": step_skip,
                "diffusivity_cm2_s": r.diffusivity,
                "diffusivity_std_cm2_s": r.diffusivity_std,
                "diffusivity_relative_standard_deviation": r.diffusivity_relative_std,
                "diffusivity_components_cm2_s": r.diffusivity_components,
                "n_jump": r.n_jump,
                "n_jump_component": r.n_jump_component,
                "fit_start_time_ps": r.fit_start_time_ps,
                "fit_end_time_ps": r.fit_end_time_ps,
                "msd_start_A2": r.msd_start,
                "msd_middle_A2": r.msd_middle,
                "msd_last_A2": r.msd_last,
                "plot_path": str(plot_path) if plot_path else "",
            }
            rows.append(row)
        pd.DataFrame(rows).to_csv(summary_csv_path, index=False)

    return DiffusivityResult(
        temperature=temperature,
        timestep_fs=timestep_fs,
        step_skip=step_skip,
        primary_species=species,
        species_results=species_results,
        msd_csv_path=msd_csv_path,
        summary_csv_path=summary_csv_path,
        plot_path=plot_path,
    )
