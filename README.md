# OxideSSE

**OxideSSE** is a Python toolkit for any type of oxide solid-state electrolyte (SSE) research. It provides reusable workflows for structure relaxation, stability screening, LAMMPS trajectory analysis, diffusivity calculation, and Arrhenius extrapolation. This repository will provide the example with a focus on Li7La3Zr2O12(LLZO)-derived garnet-type oxide materials. 

```python
import oxidesse as sse
```

## ✨ Main features

- Geometry optimization using ASE-compatible machine learning interatomic potential (MLIP) calculators.
- Binary-oxide-referenced formation energy with MLIP-recalculated reference energies.
- Energy above hull using Materials Project entries with MLIP-recalculated entry energies.
- Li-ion MSD, diffusivity, and ionic conductivity analysis from LAMMPS simulation outputs.
- Arrhenius fitting, 298 K diffusivity extrapolation, and ionic conductivity estimation.

![ToC](assets/ToC.png)

---

## 📦 Installation

### Option 1. Clone and install in editable mode

Use this option if you want to use the functions and modify the source code.

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
pip install -e .
```

Check the installation:

```bash
python -c "import oxidesse as sse; print(sse.__version__)"
```

### Option 2. Reproducible conda environment

For better reproducibility, you can create a conda environment from the provided `environment.yml` file. The `environment.yml` file records the versions used during package testing.

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
conda env create -f environment.yml
conda activate oxidesse
pip install -e .
```

---

## 🧠 SevenNet installation and requirement

OxideSSE utilizes machine-learning interatomic potentials for **geometry optimization** and **stability analysis**. The default workflow is based on **SevenNet**, but SevenNet is **not installed automatically** by OxideSSE because GPU, CUDA, PyTorch, and model compatibility can depend strongly on the user's computing environment.

> Install SevenNet separately before using SevenNet-based calculators by referring official documents: 
> - https://github.com/MDIL-SNU/SevenNet
> - https://sevennet.readthedocs.io/en/latest/

We have tested SevenNet and PyTorch version as like below with on GPU of CUDA version=13.2:
| Component | Tested version |
|---|---:|
| sevenn | 0.12.1 |
| torch | 2.11.0 |
| torch-geometric | 2.7.0 |
| torchaudio | 2.11.0 |
| torchvision | 0.26.0 |

---

Supported calculator names : `7net-0`(default), `7net-l3i5`, `7net-omat`, `7net-mf-ompa`, `7net-omni`

**Tip:** `7net-mf-ompa` and `7net-omni` are multi-modal MLIP. They need specific modality.

- Advanced users can also pass any already-created ASE calculator object directly to OxideSSE functions.

---

## 🔑 Materials Project API key

Stability calculations (formation energy and energy-above-hull calculations) require access to the Materials Project API.

Set it as an environment variable:

```bash
export MP_API_KEY="YOUR_MP_API_KEY"
```

You can also pass the API key directly in function:

```python
form = sse.compute_binary_oxide_formation_energy(
    structure="input.cif",
    api_key="YOUR_MP_API_KEY",
)

hull = sse.compute_energy_above_hull(
    structure="input.cif",
    api_key="YOUR_MP_API_KEY",
)
```

---

# 🚀 Usage with Python API

## 1. ⚙️ `sse.optimize_structure()`

Optimize a CIF or POSCAR/CONTCAR-style structure file using an ASE-compatible calculator.

### Minimal example

```python
result = sse.optimize_structure(
    structure="./examples/structures/Tetragonal_LLCO.cif",
    output_dir="./optimized"
)
```

### Variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF, POSCAR, and CONTCAR-style files are supported. |
| `output_dir` | Directory where the optimized CIF and optimization log are saved. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `output_filename` | Name of the optimized structure file. If `None`, OxideSSE uses `opt_<input>.cif`. |
| `fmax` | Force convergence criterion in eV/Å. |
| `max_steps` | Maximum number of optimizer steps. |
| `cell_relax` | If `True`, relaxes both cell and atomic positions. If `False`, relaxes atomic positions only. |
| `device` | Device option passed to the calculator builder, for example `"auto"`, `"cpu"`, or `"cuda"`. |

### Output files

- Optimized CIF file, for example `optimized/Tetragonal_LLCO.cif`
- Optimization log file, for example `optimized/Tetragonal_LLCO.opt.log`

### Result attributes

```python
result.converged
result.energy
result.energy_per_atom
result.max_force
result.output_path
result.log_path
result.message
```

---

## 2. 🧱 `sse.compute_binary_oxide_formation_energy()`

Compute the **formation energy** of an oxide structure using binary oxide references. OxideSSE automatically balances the target composition against binary oxide references, fetches representative structures from Materials Project, recalculates their energies using the selected MLIP calculator, and computes the formation energy. For the details, refer this paper. [Ong et al., Nat. Commun., 9, 3800 (2018)](https://www.nature.com/articles/s41467-018-06322-x)

Example composition balance:

```text
Li7La3Zr2O12 = 3.5 Li2O + 1.5 La2O3 + 2 ZrO2
```

### Minimal example

```python
result = sse.compute_binary_oxide_formation_energy(
    structure="./examples/structures/Tetragonal_LLCO.cif",
)
```

### Variables

| Variable | Description |
|---|---|
| `structure` | Input oxide structure. CIF and POSCAR/CONTCAR-style files are supported. The structure must contain oxygen. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY` from the environment. |
| `output_csv` | CSV path where the result is saved. Set `None` to skip CSV output. |
| `custom_reference_formulas` | Optional mapping for unsupported elements, for example `{"Mg": "MgO"}`. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Output files

- CSV file contains the target formula, target energy, formation energy, formation energy per atom, and the binary oxide references used in the calculation.

### Result attributes

```python
result.formula
result.target_energy
result.target_energy_per_atom
result.formation_energy
result.formation_energy_per_atom
result.csv_path
result.references
```

---

## 3. 📉 `sse.compute_energy_above_hull()`

Compute **energy above hull** using Materials Project entries as phase diagram references. The MP energies are re-caclulated using MLIP for the consistent fidelity.

### Minimal example

```python
result = sse.compute_energy_above_hull(
    structure="./examples/structures/Tetragonal_LLCO.cif",
)
```

### Variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF and POSCAR/CONTCAR-style files are supported. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY`. |
| `thermo_type` | MP thermo entry set. Supported values include `"GGA_GGA+U"`, `"R2SCAN"`, and `"GGA_GGA+U_R2SCAN"`. Default: `"GGA_GGA+U"`. |
| `output_csv` | CSV file for the hull result. |
| `cache_dir` | Directory where MLIP-recalculated MP entry energies are cached. |
| `use_cache` | If `True`, reuse cached MLIP energies for previously evaluated MP entries for saving the calculation time. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Output files

- CSV file includes the formula, chemical system, thermo type, target MLIP energy, energy above hull, number of MP entries used, failed entries, and cache path.

### Result attributes

```python
result.formula
result.chemical_system
result.thermo_type
result.target_energy
result.target_energy_per_atom
result.energy_above_hull
result.n_entries
result.failed_entries
result.cache_path
result.csv_path
```

---

## 4. 🚶 `sse.compute_diffusivity_from_lammps()`

Analyze **LAMMPS outputs**, convert the LAMMPS trajectory into structures, compute MSD curves, estimate diffusivity, calculate primary-species ionic conductivity using the LAMMPS data-file structure, perform error analysis, and save plots and CSV files.

### Minimal example

```python
result = sse.compute_diffusivity_from_lammps(
    simulation_dir="./examples/lammps_run",
    dump_file="dump.traj",
    data_file="structure.data",
    timestep_fs=1.0,
    step_skip=500,
    temperature=1000,
    output_dir="./diffusivity_results"
)
```

### Variables

| Variable | Description |
|---|---|
| `simulation_dir` | Directory containing the LAMMPS simulation files. |
| `dump_file` | LAMMPS dump trajectory file. Absolute paths and paths relative to `simulation_dir` are supported. |
| `data_file` | LAMMPS data file. Absolute paths and paths relative to `simulation_dir` are supported. This file is also used to obtain the cell volume and ion concentration for conductivity conversion. |
| `timestep_fs` | LAMMPS MD timestep in femtoseconds. For example, `1.0` means one LAMMPS step is 1 fs. |
| `step_skip` | Dump interval in MD steps. If frames are saved every 500 steps and `timestep_fs=1.0`, the frame spacing is 500 fs = 0.5 ps. |
| `temperature` | MD temperature in K. |
| `save_poscars` | If `True`, save converted trajectory frames as POSCAR files. Default = `False` |
| `save_plot` | If `True`, save the MSD plot. Default = `True`|
| `plot_other_species` | If `True`, include non-primary species in the MSD plot. Default = `True`|
| `species` | Primary mobile species. Default: `"Li"`. |
| `output_dir` | Directory for MSD plot and CSV files. |
| `output_csv` | Summary CSV path. |
| `msd_csv` | Time-resolved MSD CSV path. |
| `oxidized_species` | Optional oxidation-state mapping for conductivity conversion, for example `{"Li": "Li+"}`. If omitted and `species="Li"`, Li is treated as `Li+`. |

### Output files

- MSD plot, for example `diffusivity_results/msd_<species>_<temperature>K.png`
- Diffusivity and conductivity analysis results, for example `diffusivity_results/diffusivity_summary.csv`
- MSD raw data in each species, for example `diffusivity_results/msd_by_species.csv`
- Optional converted POSCAR files if `save_poscars=True`, for example `examples/lammps_run/traj_to_POSCARs`

### Result attributes

```python
result.temperature
result.timestep_fs
result.step_skip
result.primary_species
result.plot_path
result.summary_csv_path
result.msd_csv_path
result.conductivity
result.conductivity_std
result.conductivity_source_file
result.species_results
```

---

## 5. 🔥 `sse.plot_arrhenius()`

Fit Arrhenius behavior from temperature-dependent diffusivity data, save an Arrhenius plot, extrapolate diffusivity to 298 K, and estimate 298 K ionic conductivity.

### Minimal example

```python
data = {
    1000: {'diffusivity': 7.70E-06, 'std': 4.90E-07},
    900: {'diffusivity': 5.19E-06, 'std': 3.62E-07},
    800: {'diffusivity': 3.05E-06, 'std': 2.35E-07},
    700: {'diffusivity': 1.70E-06, 'std': 1.55E-07},
    600: {'diffusivity': 1.17E-06, 'std': 1.20E-07},
}

result = sse.plot_arrhenius(
    data=data,
    name="t-LLCO",
    structure="./examples/structures/Tetragonal_LLCO.cif",
    output_dir="./arrhenius_results"
)
```

`structure` may also be a LAMMPS data file, for example `structure="./examples/lammps_run/system.data"`.

### Variables

| Variable | Description |
|---|---|
| `data` | Dictionary containing temperature, diffusivity, and standard deviation. Required format: `{T: {"diffusivity": D, "std": std}}`. |
| `name` | Label for the measured diffusivity data. |
| `structure` | Structure file used to compute 298 K ionic conductivity via the Nernst-Einstein relation. Supported formats: CIF, POSCAR/CONTCAR, and LAMMPS data file. |
| `output_dir` | Directory where the plot and CSV are saved. |
| `extrapolate_298K` | If `True`, plot the fitted line to 298 K and report 298 K diffusivity. |
| `output_filename` | Arrhenius plot filename. |
| `output_csv` | Summary CSV filename. Set `None` to skip CSV output. |
| `specie` | Oxidized species string, usually `"Li+"`. |
| `conductivity_factor_298K` | Optional precomputed conversion factor from diffusivity to conductivity. If provided, `structure` is not required for conductivity. |

### Output files

- Arrhenius plot, for example `arrhenius_results/arrhenius_plot.png`
- Summary CSV, for example `arrhenius_results/arrhenius_summary.csv`

### Result attributes

```python
result.activation_energy_eV
result.activation_energy_error_eV
result.diffusivity_298K_cm2_s
result.diffusivity_298K_std_cm2_s
result.diffusivity_298K_range_cm2_s
result.conductivity_298K_mS_cm
result.conductivity_298K_range_mS_cm
result.plot_path
result.csv_path
```

---

# 🖥️ Usage with command-line interface (CLI)

OxideSSE provides an `argparse`-based command-line interface.

```bash
oxidesse --help
```

## 1. ⚙️ Geometry optimization

```bash
oxidesse optimize-structure \
  --structure input.cif
```

Options:

- `--structure`: input CIF/POSCAR path
- `--calculator`: calculator name, default `7net-0`
- `--output_dir`: output directory
- `--output_filename`: optimized CIF filename
- `--fmax`: force convergence criterion
- `--max_steps`: maximum optimization steps
- `--no_cell_relax`: disable cell relaxation and optimize atomic positions only
- `--device`: calculator device option

## 2. 🧱 Formation energy

```bash
oxidesse formation-energy \
  --structure input.cif
```

Options:

- `--structure`: input oxide structure
- `--calculator`: calculator name, default `7net-0`
- `--api_key`: Materials Project API key, optional if `MP_API_KEY` is set
- `--output_csv`: CSV output path
- `--device`: calculator device option
- `--no_append`: overwrite CSV instead of appending to an existing file

## 3. 📉 Energy above hull

```bash
oxidesse energy-above-hull \
  --structure input.cif
```

Options:

- `--structure`: input structure
- `--calculator`: calculator name, default `7net-0`
- `--api_key`: Materials Project API key, optional if `MP_API_KEY` is set
- `--thermo_type`: `GGA_GGA+U` (default), `GGAU`, `R2SCAN`, `GGA_GGA+U_R2SCAN`, or `GGAU_R2SCAN`
- `--output_csv`: CSV output path
- `--cache_dir`: cache directory
- `--no_use_cache`: disable reuse of cached MLIP entry energies
- `--no_append`: overwrite CSV instead of appending to an existing file

## 4. 🚶 Diffusivity from LAMMPS

```bash
oxidesse diffusivity \
  --simulation_dir ./examples/lammps_run \
  --dump_file dump.traj \
  --data_file structure.data \
  --timestep_fs 1.0 \
  --step_skip 500 \
  --temperature 1000
```

Options:

- `--simulation_dir`: LAMMPS simulation directory
- `--dump_file`: LAMMPS dump file
- `--data_file`: LAMMPS data file
- `--timestep_fs`: MD timestep in fs
- `--step_skip`: dump interval in MD steps
- `--temperature`: temperature in K
- `--species`: primary mobile species, default `Li`
- `--output_dir`: output directory
- `--save_poscars`: save converted POSCAR frames
- `--no_save_plot`: disable MSD plot saving
- `--no_plot_other_species`: plot only the primary species
- `--oxidized_specie`: oxidized primary species for conductivity conversion, for example `Li+`

## 5. 🔥 Arrhenius fitting

For command-line use, pass data as a Python-style dictionary string:

```bash
oxidesse arrhenius \
  --data '{1000: {"diffusivity": 7.70E-06, "std": 4.90E-07}, 900: {"diffusivity": 5.19E-06, "std": 3.62E-07}, 800: {"diffusivity": 3.05E-06, "std": 2.35E-07}, 700: {"diffusivity": 1.70E-06, "std": 1.55E-07}, 600: {"diffusivity": 1.17E-06, "std": 1.20E-07}}' \
  --name t-LLCO \
  --structure input.cif
```

Options:

- `--data`: temperature-diffusivity dictionary
- `--name`: dataset label, title of arrhenius plot
- `--output_dir`: output directory
- `--output_filename`: plot filename
- `--output_csv`: CSV filename
- `--extrapolate_298K`: whether to extrapolate to 298 K
- `--structure`: structure used for conductivity conversion. Supported formats: CIF, POSCAR/CONTCAR, and LAMMPS data file
- `--specie`: oxidized mobile species, for example `Li+`
- `--conductivity_factor_298K`: optional precomputed conductivity conversion factor

---

# 📁 Repository layout

```text
OxideSSE/
├── src/oxidesse/
│   ├── _aimd.py
│   ├── arrhenius.py
│   ├── binary_oxide.py
│   ├── calculators.py
│   ├── diffusion.py
│   ├── formation_energy.py
│   ├── hull.py
│   ├── io.py
│   └── optimization.py
├── examples/
├── assets/
├── environment.yml
├── pyproject.toml
└── README.md
```

---

# 📚 Citation

If you use OxideSSE in your research, please cite:

```bibtex
@misc{jeon2026unlocking,
  title = {Harnessing Universal Machine Learning Interatomic Potentials for Unlocking the Li7M3X2O12 Garnet Electrolyte Landscape},
  author = {Jeon, Mingyu and Lee, Jae-Kwan and Artrith, Nongnuch and Urban, Alexander and Lee, Byungju and Kim, Jieun and Kim, Donghun and Lee, Jung-Hoon},
  year = {2026},
  archivePrefix = {arXiv},
  eprint = {},
  url = {}
}
```

---

# 📄 License

OxideSSE is distributed under the MIT License. See the [`LICENSE`](LICENSE) file for details.
