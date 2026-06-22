"""Arrhenius analysis for diffusivity and room-temperature extrapolation."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

K_B_EV_PER_K = 8.617333262145e-5
def _default_color():
    import matplotlib.pyplot as plt
    return tuple(plt.cm.Blues(np.linspace(0.95, 0.3, 4))[0])

DEFAULT_COLOR = None


@dataclass
class ArrheniusResult:
    temperatures_K: list[float]
    diffusivities_cm2_s: list[float]
    std_cm2_s: list[float]
    slope: float
    intercept: float
    activation_energy_eV: float
    activation_energy_error_eV: float | None
    log10_diffusivity_298K_std: float | None
    diffusivity_298K_std_cm2_s: float | None
    diffusivity_298K_cm2_s: float | None
    diffusivity_298K_range_cm2_s: tuple[float, float] | None
    conductivity_298K_mS_cm: float | None
    conductivity_298K_range_mS_cm: tuple[float, float] | None
    plot_path: Path | None = None
    csv_path: Path | None = None


def _normalize_data(data: Mapping[float, Mapping[str, float]]):
    temps = sorted([float(t) for t in data.keys()], reverse=True)
    diffusivities = []
    stds = []
    for t in temps:
        item = data[t] if t in data else data[int(t)]
        diff = item.get("diffusivity", item.get("diff", item.get("D")))
        std = item.get("std", item.get("diffusivity_std", item.get("error", 0.0)))
        if diff is None:
            raise ValueError("Each temperature entry must contain 'diffusivity'.")
        if diff <= 0:
            raise ValueError("Diffusivity values must be positive for Arrhenius fitting.")
        diffusivities.append(float(diff))
        stds.append(float(std or 0.0))
    if len(temps) < 2:
        raise ValueError("At least two temperature points are required for Arrhenius fitting.")
    return np.asarray(temps), np.asarray(diffusivities), np.asarray(stds)


def get_conductivity_conversion_factor(structure, specie: str = "Li+", temperature: float = 298.0) -> float:
    """Convert diffusivity in cm^2/s to conductivity in mS/cm.

    The conversion uses the Nernst-Einstein equation and therefore requires a
    structure volume, species concentration, and oxidation-decorated specie such
    as ``"Li+"``.
    """
    from pymatgen.core.periodic_table import Specie

    if not hasattr(structure, "composition"):
        from .io import read_structure

        structure = read_structure(structure)

    df_sp = specie if isinstance(specie, Specie) else Specie.from_str(specie)
    z = df_sp.oxi_state
    n = structure.composition.get(str(df_sp.element), 0)
    if n == 0:
        n = structure.composition.get(df_sp, 0)
    if n == 0:
        raise ValueError(f"No {specie} found in structure composition: {structure.composition}")

    volume_cm3 = structure.volume * 1e-24
    avogadro = 6.022140857e23
    elementary_charge = 1.6021766208e-19
    gas_constant = 8.3144598
    return 1000 * n / (volume_cm3 * avogadro) * z**2 * (avogadro * elementary_charge) ** 2 / (gas_constant * temperature)


def fit_arrhenius(
    data: Mapping[float, Mapping[str, float]],
    extrapolate_298K: bool = True,
    structure: Any | None = None,
    specie: str = "Li+",
    conductivity_factor_298K: float | None = None,
) -> ArrheniusResult:
    """Fit log10(D) versus 1000/T and optionally extrapolate to 298 K."""
    temps, diffusivities, stds = _normalize_data(data)
    x = 1000.0 / temps
    y = np.log10(diffusivities)
    yerr = np.zeros_like(stds, dtype=float)
    mask = (stds > 0) & (diffusivities > 0)
    yerr[mask] = stds[mask] / (diffusivities[mask] * np.log(10))

    def _linear_model(x_values, slope_value, intercept_value):
        return slope_value * x_values + intercept_value

    # Use scipy.optimize.curve_fit to stay consistent with the original
    # notebook workflow. For a linear Arrhenius model, curve_fit and
    # np.polyfit solve essentially the same least-squares problem, but
    # curve_fit makes the covariance handling explicit and matches the
    # research notebook implementation.
    if np.all(yerr > 0):
        popt, cov = curve_fit(
            _linear_model,
            x,
            y,
            sigma=yerr,
            absolute_sigma=True,
            maxfev=10000,
        )
    else:
        popt, cov = curve_fit(
            _linear_model,
            x,
            y,
            maxfev=10000,
        )
    slope, intercept = float(popt[0]), float(popt[1])
    slope_err = float(np.sqrt(cov[0, 0])) if cov is not None and cov.size else None
    ea = -slope * K_B_EV_PER_K * np.log(10) * 1000.0
    ea_err = None if slope_err is None else slope_err * K_B_EV_PER_K * np.log(10) * 1000.0

    logd298_std = None
    d298_std = None
    d298 = None
    d298_range = None
    cond298 = None
    cond298_range = None

    if extrapolate_298K:
        x298 = 1000.0 / 298.0
        logd298 = slope * x298 + intercept
        d298 = float(10**logd298)

        # Propagate the uncertainty of the fitted Arrhenius line to the
        # extrapolated 298 K point.  The fit is performed in log10(D) versus
        # 1000/T space, so the standard deviation is first estimated in
        # log10(D) and then converted into an asymmetric diffusivity range.
        if cov is not None and np.size(cov) == 4:
            design = np.array([x298, 1.0], dtype=float)
            variance = float(design @ cov @ design.T)
            if variance > 0:
                logd298_std = float(np.sqrt(variance))
                d298_min = float(10 ** (logd298 - logd298_std))
                d298_max = float(10 ** (logd298 + logd298_std))
                d298_range = (d298_min, d298_max)
                d298_std = float((d298_max - d298_min) / 2.0)
            else:
                d298_range = (d298, d298)
                d298_std = 0.0
        elif slope_err is not None:
            # Fallback for environments where covariance is unavailable.
            logd_min = (slope - slope_err) * x298 + intercept
            logd_max = (slope + slope_err) * x298 + intercept
            d298_range = (float(min(10**logd_min, 10**logd_max)), float(max(10**logd_min, 10**logd_max)))
            d298_std = float((d298_range[1] - d298_range[0]) / 2.0)
        else:
            d298_range = (d298, d298)
            d298_std = 0.0

        factor = conductivity_factor_298K
        if factor is None and structure is not None:
            factor = get_conductivity_conversion_factor(structure, specie=specie, temperature=298.0)
        if factor is not None:
            cond298 = float(factor * d298)
            cond298_range = (float(factor * d298_range[0]), float(factor * d298_range[1]))

    return ArrheniusResult(
        temperatures_K=list(map(float, temps)),
        diffusivities_cm2_s=list(map(float, diffusivities)),
        std_cm2_s=list(map(float, stds)),
        slope=slope,
        intercept=intercept,
        activation_energy_eV=float(ea),
        activation_energy_error_eV=ea_err,
        log10_diffusivity_298K_std=logd298_std,
        diffusivity_298K_std_cm2_s=d298_std,
        diffusivity_298K_cm2_s=d298,
        diffusivity_298K_range_cm2_s=d298_range,
        conductivity_298K_mS_cm=cond298,
        conductivity_298K_range_mS_cm=cond298_range,
    )


def plot_arrhenius(
    data: Mapping[float, Mapping[str, float]],
    output_dir: str | Path = ".",
    output_filename: str = "arrhenius_plot.png",
    output_csv: str | Path | None = "arrhenius_summary.csv",
    name: str = "structure",
    extrapolate_298K: bool = True,
    structure: Any | None = None,
    specie: str = "Li+",
    conductivity_factor_298K: float | None = None,
    color=None,
    linestyle: str = "-",
    marker: str = "D",
    dpi: int = 300,
) -> ArrheniusResult:
    """Fit Arrhenius data, save a plot, and write a summary CSV.

    ``data`` should be a dictionary such as ``{600: {"diffusivity": 1e-6,
    "std": 1e-7}, 700: {"diffusivity": 2e-6, "std": 2e-7}}``.
    """
    import matplotlib.pyplot as plt

    if color is None:
        color = _default_color()

    result = fit_arrhenius(
        data=data,
        extrapolate_298K=extrapolate_298K,
        structure=structure,
        specie=specie,
        conductivity_factor_298K=conductivity_factor_298K,
    )

    output_dir = Path(output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / output_filename

    temps = np.asarray(result.temperatures_K)
    diffs = np.asarray(result.diffusivities_cm2_s)
    stds = np.asarray(result.std_cm2_s)
    x = 1000.0 / temps

    fig, ax = plt.subplots(figsize=(5, 4))
    positive_std = np.where(stds > 0, stds, None)
    ax.errorbar(
        x,
        diffs,
        yerr=positive_std,
        color=color,
        marker=marker,
        linestyle="None",
        markerfacecolor=color,
        markeredgecolor="black",
        markersize=7,
        capsize=5,
        label=name,
    )

    x_end = 1000.0 / 298.0 if extrapolate_298K else max(x)
    x_line = np.linspace(min(x), x_end, 200)
    y_line = 10 ** (result.slope * x_line + result.intercept)
    ax.plot(x_line, y_line, linestyle=linestyle, color=color, linewidth=2)

    if extrapolate_298K and result.diffusivity_298K_cm2_s is not None:
        x298 = 1000.0 / 298.0
        y298 = result.diffusivity_298K_cm2_s
        if result.diffusivity_298K_range_cm2_s is not None:
            yerr_298 = np.array([[
                max(y298 - result.diffusivity_298K_range_cm2_s[0], 0.0),
                max(result.diffusivity_298K_range_cm2_s[1] - y298, 0.0),
            ]]).T
            ax.errorbar(
                [x298],
                [y298],
                yerr=yerr_298,
                color=color,
                marker=marker,
                linestyle="None",
                markerfacecolor=color,
                markeredgecolor="black",
                markersize=7,
                capsize=5,
                zorder=6,
            )
        else:
            ax.scatter([x298], [y298], color=color, marker=marker, edgecolor="black", zorder=5)

    ax.set_yscale("log")
    ax.set_xlabel(r"1000 / T (K$^{-1}$)")
    ax.set_ylabel(r"Diffusivity (cm$^2$/s)")
    ax.set_title(f"Arrhenius plot: {name}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=dpi)
    plt.close(fig)
    result.plot_path = plot_path

    if output_csv is not None:
        csv_path = Path(output_csv)
        if not csv_path.is_absolute():
            csv_path = output_dir / csv_path
        _write_arrhenius_csv(result, csv_path, name=name)
        result.csv_path = csv_path
    return result


def _write_arrhenius_csv(result: ArrheniusResult, output_csv: Path, name: str) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "name": name,
        "temperatures_K": ";".join(f"{t:g}" for t in result.temperatures_K),
        "diffusivities_cm2_s": ";".join(f"{d:.12g}" for d in result.diffusivities_cm2_s),
        "std_cm2_s": ";".join(f"{s:.12g}" for s in result.std_cm2_s),
        "slope": result.slope,
        "intercept": result.intercept,
        "activation_energy_eV": result.activation_energy_eV,
        "activation_energy_error_eV": result.activation_energy_error_eV,
        "diffusivity_298K_cm2_s": result.diffusivity_298K_cm2_s,
        "diffusivity_298K_std_cm2_s": result.diffusivity_298K_std_cm2_s,
        "log10_diffusivity_298K_std": result.log10_diffusivity_298K_std,
        "diffusivity_298K_min_cm2_s": None if result.diffusivity_298K_range_cm2_s is None else result.diffusivity_298K_range_cm2_s[0],
        "diffusivity_298K_max_cm2_s": None if result.diffusivity_298K_range_cm2_s is None else result.diffusivity_298K_range_cm2_s[1],
        "conductivity_298K_mS_cm": result.conductivity_298K_mS_cm,
        "conductivity_298K_min_mS_cm": None if result.conductivity_298K_range_mS_cm is None else result.conductivity_298K_range_mS_cm[0],
        "conductivity_298K_max_mS_cm": None if result.conductivity_298K_range_mS_cm is None else result.conductivity_298K_range_mS_cm[1],
        "plot_path": str(result.plot_path) if result.plot_path else "",
    }
    pd.DataFrame([row]).to_csv(output_csv, index=False)
