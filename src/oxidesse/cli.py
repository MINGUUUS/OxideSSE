"""Command-line interface for OxideSSE."""
from __future__ import annotations

import argparse
import ast
import os
from typing import Any

from .arrhenius import plot_arrhenius
from .diffusion import compute_diffusivity_from_lammps
from .formation_energy import compute_binary_oxide_formation_energy
from .hull import compute_energy_above_hull_mlp
from .optimization import optimize_structure, optimize_structures


def _none_if_string_none(value: str | None) -> str | None:
    if value is None:
        return None
    return None if str(value).strip().lower() == "none" else value


def _parse_mapping(value: str | None) -> dict[str, Any] | None:
    if value is None:
        return None
    try:
        parsed = ast.literal_eval(value)
    except Exception as exc:
        raise argparse.ArgumentTypeError(
            "Expected a Python/JSON-style dictionary string, for example "
            "'{\"Li\": \"Li+\"}'."
        ) from exc
    if not isinstance(parsed, dict):
        raise argparse.ArgumentTypeError("Expected a dictionary.")
    return parsed


def _print_path(label: str, path) -> None:
    if path:
        print(f"{label}: {path}")


def _print_diffusivity_result(result) -> None:
    print("Diffusivity analysis completed.")
    print(f"Primary species: {result.primary_species}")
    print(f"Temperature: {result.temperature:g} K")
    if result.diffusivity is not None:
        print(f"Diffusivity: {result.diffusivity:.8e} cm^2/s")
    if result.diffusivity_std is not None:
        print(f"Diffusivity std: {result.diffusivity_std:.8e} cm^2/s")
    if result.conductivity is not None:
        print(f"Conductivity: {result.conductivity:.8e} mS/cm")
    if result.conductivity_std is not None:
        print(f"Conductivity std: {result.conductivity_std:.8e} mS/cm")
    _print_path("Summary CSV", result.summary_csv_path)
    _print_path("MSD CSV", result.msd_csv_path)
    _print_path("MSD plot", result.plot_path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="oxidesse",
        description="OxideSSE command-line interface.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser(
        "optimize-structure",
        help="Optimize one CIF/POSCAR structure.",
        description="Optimize one CIF/POSCAR structure using an ASE-compatible calculator.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--structure", required=True, metavar="PATH", help="Input CIF/POSCAR/CONTCAR structure file. Type: path string.")
    p.add_argument("--calculator", default="7net-0", metavar="NAME", help="ASE calculator name or supported SevenNet model name. Type: string.")
    p.add_argument("--output_dir", "--output-dir", dest="output_dir", default=".", metavar="DIR", help="Directory where optimized structure and log files are saved. Type: path string.")
    p.add_argument("--output_filename", "--output-filename", dest="output_filename", default=None, metavar="FILENAME", help="Output optimized structure filename. Type: string or None.")
    p.add_argument("--fmax", type=float, default=0.01, metavar="FLOAT", help="Force convergence criterion in eV/Angstrom. Type: float.")
    p.add_argument("--max_steps", "--max-steps", dest="max_steps", type=int, default=500, metavar="INT", help="Maximum number of optimizer steps. Type: int.")
    p.add_argument("--no_cell_relax", "--no-cell-relax", dest="cell_relax", action="store_false", default=True, help="Disable cell relaxation and optimize atomic positions only. Type: boolean flag.")
    p.add_argument("--device", default="auto", metavar="DEVICE", help="Calculator device, for example auto, cpu, cuda, or cuda:0. Type: string.")

    p = subparsers.add_parser(
        "optimize-structures",
        help="Batch optimize structures in a directory.",
        description="Batch optimize structures matching a file pattern in an input directory.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--input_dir", "--input-dir", dest="input_dir", required=True, metavar="DIR", help="Directory containing input structure files. Type: path string.")
    p.add_argument("--pattern", default="*.cif", metavar="GLOB", help="File pattern used to select structures, for example *.cif. Type: string.")
    p.add_argument("--output_dir", "--output-dir", dest="output_dir", default="optimized", metavar="DIR", help="Directory where optimized structures are saved. Type: path string.")
    p.add_argument("--calculator", default="7net-0", metavar="NAME", help="ASE calculator name or supported SevenNet model name. Type: string.")
    p.add_argument("--fmax", type=float, default=0.01, metavar="FLOAT", help="Force convergence criterion in eV/Angstrom. Type: float.")
    p.add_argument("--max_steps", "--max-steps", dest="max_steps", type=int, default=500, metavar="INT", help="Maximum number of optimizer steps. Type: int.")
    p.add_argument("--no_cell_relax", "--no-cell-relax", dest="cell_relax", action="store_false", default=True, help="Disable cell relaxation and optimize atomic positions only. Type: boolean flag.")
    p.add_argument("--device", default="auto", metavar="DEVICE", help="Calculator device, for example auto, cpu, cuda, or cuda:0. Type: string.")

    p = subparsers.add_parser(
        "formation-energy",
        help="Compute binary-oxide-referenced formation energy.",
        description="Compute binary-oxide-referenced formation energy using MLP-recalculated target and reference energies.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--structure", required=True, metavar="PATH", help="Input oxide structure file. Type: path string.")
    p.add_argument("--calculator", default="7net-0", metavar="NAME", help="ASE calculator name or supported SevenNet model name. Type: string.")
    p.add_argument("--api_key", "--api-key", dest="api_key", default=None, metavar="KEY", help="Materials Project API key. If omitted, MP_API_KEY is used. Type: string or None.")
    p.add_argument("--output_csv", "--output-csv", dest="output_csv", default=None, metavar="PATH", help="CSV file path for saving the result. Use None to skip CSV output. Type: path string or None.")
    p.add_argument("--device", default="auto", metavar="DEVICE", help="Calculator device, for example auto, cpu, cuda, or cuda:0. Type: string.")
    p.add_argument("--no_append", "--no-append", dest="append", action="store_false", default=True, help="Overwrite CSV instead of appending to an existing file. Type: boolean flag.")

    p = subparsers.add_parser(
        "energy-above-hull",
        help="Compute MLP-recalculated energy above hull.",
        description="Compute energy above hull using Materials Project entries with MLP-recalculated energies.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--structure", required=True, metavar="PATH", help="Input structure file. Type: path string.")
    p.add_argument("--calculator", default="7net-0", metavar="NAME", help="ASE calculator name or supported SevenNet model name. Type: string.")
    p.add_argument("--api_key", "--api-key", dest="api_key", default=None, metavar="KEY", help="Materials Project API key. If omitted, MP_API_KEY is used. Type: string or None.")
    p.add_argument("--thermo_type", "--thermo-type", dest="thermo_type", default="GGA_GGA+U", metavar="TYPE", help="Materials Project thermo entry set. Supported: GGA_GGA+U, R2SCAN, or GGA_GGA+U_R2SCAN. Default: GGA_GGA+U. Type: string.")
    p.add_argument("--output_csv", "--output-csv", dest="output_csv", default=None, metavar="PATH", help="CSV file path for saving the result. Use None to skip CSV output. Type: path string or None.")
    p.add_argument("--cache_dir", "--cache-dir", dest="cache_dir", default=".oxidesse_cache", metavar="DIR", help="Directory for cached MLP-recalculated MP entry energies. Use None to disable caching. Type: path string or None.")
    p.add_argument("--no_use_cache", "--no-use-cache", dest="use_cache", action="store_false", default=True, help="Disable reuse of cached MLP energies. Type: boolean flag.")
    p.add_argument("--device", default="auto", metavar="DEVICE", help="Calculator device, for example auto, cpu, cuda, or cuda:0. Type: string.")
    p.add_argument("--no_append", "--no-append", dest="append", action="store_false", default=True, help="Overwrite CSV instead of appending to an existing file. Type: boolean flag.")

    p = subparsers.add_parser(
        "diffusivity",
        help="Compute MSD, diffusivity, and conductivity from a LAMMPS simulation.",
        description="Compute MSD, diffusivity, error estimates, and primary-species ionic conductivity from one LAMMPS dump/data pair.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--simulation_dir", "--simulation-dir", dest="simulation_dir", required=True, metavar="DIR", help="Directory containing the LAMMPS dump and data files. Type: path string.")
    p.add_argument("--dump_file", "--dump-file", dest="dump_file", required=True, metavar="PATH", help="LAMMPS dump trajectory file. Absolute paths or paths relative to simulation_dir are supported. Type: path string.")
    p.add_argument("--data_file", "--data-file", dest="data_file", required=True, metavar="PATH", help="LAMMPS data file used for topology, cell volume, species, and conductivity conversion. Absolute paths or paths relative to simulation_dir are supported. Type: path string.")
    p.add_argument("--timestep_fs", "--timestep-fs", dest="timestep_fs", required=True, type=float, metavar="FLOAT", help="LAMMPS MD timestep in femtoseconds. Type: float.")
    p.add_argument("--step_skip", "--step-skip", dest="step_skip", required=True, type=int, metavar="INT", help="Dump interval in MD steps. Type: int.")
    p.add_argument("--temperature", required=True, type=float, metavar="FLOAT", help="Simulation temperature in K, used for diffusivity analysis and conductivity conversion. Type: float.")
    p.add_argument("--species", default="Li", metavar="ELEMENT", help="Primary mobile species for diffusivity and conductivity output. Type: element symbol string.")
    p.add_argument("--output_dir", "--output-dir", dest="output_dir", default=".", metavar="DIR", help="Directory where MSD plot and CSV files are saved. Type: path string.")
    p.add_argument("--save_poscars", "--save-poscars", dest="save_poscars", action="store_true", default=False, help="Save converted trajectory frames as POSCAR files. Type: boolean flag.")
    p.add_argument("--no_save_plot", "--no-save-plot", dest="save_plot", action="store_false", default=True, help="Disable saving the MSD plot. Type: boolean flag.")
    p.add_argument("--no_plot_other_species", "--no-plot-other-species", dest="plot_other_species", action="store_false", default=True, help="Plot only the primary species and omit other species from the MSD plot. Type: boolean flag.")
    p.add_argument("--output_csv", "--output-csv", dest="output_csv", default="diffusivity_summary.csv", metavar="PATH", help="Summary CSV filename or path. Use None to skip CSV output. Type: path string or None.")
    p.add_argument("--msd_csv", "--msd-csv", dest="msd_csv", default="msd_by_species.csv", metavar="PATH", help="Time-resolved MSD CSV filename or path. Use None to skip MSD CSV output. Type: path string or None.")
    p.add_argument("--oxidized_specie", "--oxidized-specie", dest="oxidized_specie", default=None, metavar="SPECIE", help="Oxidized primary species used for conductivity conversion, for example Li+. If omitted and species=Li, Li+ is assumed. Type: string or None.")
    p.add_argument("--oxidized_species", "--oxidized-species", dest="oxidized_species", default=None, type=_parse_mapping, metavar="DICT", help="Optional oxidation-state mapping, for example '{\"Li\": \"Li+\"}'. Type: dictionary string or None.")

    p = subparsers.add_parser(
        "arrhenius",
        help="Fit Arrhenius diffusivity data.",
        description="Fit log10(D) versus 1000/T, extrapolate diffusivity to 298 K, and optionally compute 298 K conductivity.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    p.add_argument("--data", required=True, metavar="DICT", help="Temperature-diffusivity dictionary, e.g. '{1000: {\"diffusivity\": 1e-6, \"std\": 1e-7}}'. Type: dictionary string.")
    p.add_argument("--output_dir", "--output-dir", dest="output_dir", default=".", metavar="DIR", help="Directory where Arrhenius plot and CSV are saved. Type: path string.")
    p.add_argument("--output_filename", "--output-filename", dest="output_filename", default="arrhenius_plot.png", metavar="FILENAME", help="Arrhenius plot filename. Type: string.")
    p.add_argument("--output_csv", "--output-csv", dest="output_csv", default="arrhenius_summary.csv", metavar="PATH", help="Summary CSV filename or path. Use None to skip CSV output. Type: path string or None.")
    p.add_argument("--name", default="structure", metavar="LABEL", help="Dataset label used in the plot and CSV. Type: string.")
    p.add_argument("--no_extrapolate_298K", "--no-extrapolate-298K", dest="extrapolate_298K", action="store_false", default=True, help="Disable extrapolation of the fitted Arrhenius line to 298 K. Type: boolean flag.")
    p.add_argument("--structure", default=None, metavar="PATH", help="Structure file used for 298 K conductivity conversion. Supported: CIF, POSCAR/CONTCAR, or LAMMPS data file. Type: path string or None.")
    p.add_argument("--specie", default="Li+", metavar="SPECIE", help="Oxidized mobile species for conductivity conversion, for example Li+. Type: string.")
    p.add_argument("--conductivity_factor_298K", "--conductivity-factor-298K", dest="conductivity_factor_298K", type=float, default=None, metavar="FLOAT", help="Precomputed diffusivity-to-conductivity conversion factor at 298 K. Type: float or None.")

    return parser


def main(argv: list[str] | None = None) -> None:
    os.environ.setdefault("PAGER", "cat")

    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "optimize-structure":
        result = optimize_structure(
            structure=args.structure,
            calculator=args.calculator,
            output_dir=args.output_dir,
            output_filename=args.output_filename,
            fmax=args.fmax,
            max_steps=args.max_steps,
            cell_relax=args.cell_relax,
            device=args.device,
        )
        print("Geometry optimization completed." if result.converged else "Geometry optimization finished without convergence.")
        if result.energy is not None:
            print(f"Energy: {result.energy:.8f} eV")
        _print_path("Output structure", result.output_path)
        _print_path("Log file", result.log_path)
        return

    if args.command == "optimize-structures":
        results = optimize_structures(
            input_dir=args.input_dir,
            pattern=args.pattern,
            output_dir=args.output_dir,
            calculator=args.calculator,
            fmax=args.fmax,
            max_steps=args.max_steps,
            cell_relax=args.cell_relax,
            device=args.device,
        )
        n_ok = sum(1 for r in results if r.converged)
        print(f"Batch optimization completed: {n_ok}/{len(results)} converged.")
        return

    if args.command == "formation-energy":
        result = compute_binary_oxide_formation_energy(
            structure=args.structure,
            calculator=args.calculator,
            api_key=args.api_key,
            output_csv=_none_if_string_none(args.output_csv),
            device=args.device,
            append=args.append,
        )
        print("Formation-energy calculation completed.")
        print(f"Formula: {result.formula}")
        print(f"Formation energy: {result.formation_energy:.8f} eV")
        print(f"Formation energy per atom: {result.formation_energy_per_atom:.8f} eV/atom")
        _print_path("CSV", result.csv_path)
        return

    if args.command == "energy-above-hull":
        result = compute_energy_above_hull_mlp(
            structure=args.structure,
            calculator=args.calculator,
            api_key=args.api_key,
            thermo_type=args.thermo_type,
            output_csv=_none_if_string_none(args.output_csv),
            cache_dir=_none_if_string_none(args.cache_dir),
            use_cache=args.use_cache,
            device=args.device,
            append=args.append,
        )
        print("Energy-above-hull calculation completed.")
        print(f"Formula: {result.formula}")
        print(f"Chemical system: {result.chemical_system}")
        print(f"Energy above hull: {result.energy_above_hull:.8f} eV/atom")
        print(f"MP entries recalculated: {result.n_entries}")
        _print_path("CSV", result.csv_path)
        _print_path("Cache", result.cache_path)
        return

    if args.command == "diffusivity":
        oxidized_species = args.oxidized_species or None
        if args.oxidized_specie is not None:
            oxidized_species = dict(oxidized_species or {})
            oxidized_species[args.species] = args.oxidized_specie
        result = compute_diffusivity_from_lammps(
            simulation_dir=args.simulation_dir,
            dump_file=args.dump_file,
            data_file=args.data_file,
            timestep_fs=args.timestep_fs,
            step_skip=args.step_skip,
            temperature=args.temperature,
            species=args.species,
            output_dir=args.output_dir,
            save_poscars=args.save_poscars,
            save_plot=args.save_plot,
            output_csv=_none_if_string_none(args.output_csv),
            msd_csv=_none_if_string_none(args.msd_csv),
            plot_other_species=args.plot_other_species,
            oxidized_species=oxidized_species,
        )
        _print_diffusivity_result(result)
        return

    if args.command == "arrhenius":
        data = ast.literal_eval(args.data)
        result = plot_arrhenius(
            data=data,
            output_dir=args.output_dir,
            output_filename=args.output_filename,
            output_csv=_none_if_string_none(args.output_csv),
            name=args.name,
            extrapolate_298K=args.extrapolate_298K,
            structure=args.structure,
            specie=args.specie,
            conductivity_factor_298K=args.conductivity_factor_298K,
        )
        print("Arrhenius fitting completed.")
        print(f"Activation energy: {result.activation_energy_eV:.8f} eV")
        if result.diffusivity_298K_cm2_s is not None:
            print(f"Diffusivity at 298 K: {result.diffusivity_298K_cm2_s:.8e} cm^2/s")
        if result.conductivity_298K_mS_cm is not None:
            print(f"Conductivity at 298 K: {result.conductivity_298K_mS_cm:.8e} mS/cm")
        _print_path("Plot", result.plot_path)
        _print_path("CSV", result.csv_path)
        return


if __name__ == "__main__":
    main()
