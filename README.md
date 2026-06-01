# 🧪 OxideSSE

**OxideSSE** is a Python toolkit for oxide solid-state electrolyte research, with a focus on LLZO-derived garnet-type oxide materials. It provides reusable workflows for structure relaxation, stability screening, LAMMPS trajectory analysis, diffusivity calculation, and Arrhenius extrapolation.

```python
import oxidesse as sse
```

## ✨ Main features

- ⚙️ Geometry optimization using ASE-compatible MLP calculators.
- 🧱 Binary-oxide-referenced formation energy with MLP-recalculated reference energies.
- 📉 Energy above hull using Materials Project entries with MLP-recalculated entry energies.
- 🚶 Li-ion MSD and diffusivity analysis from a single LAMMPS simulation folder.
- 🔥 Arrhenius fitting, 298 K diffusivity extrapolation, and ionic conductivity estimation.

---

## 📦 Installation

### Option 1. Install directly from GitHub

```bash
pip install git+https://github.com/MINGUUUS/OxideSSE.git
```

To install a specific release tag:

```bash
pip install git+https://github.com/MINGUUUS/OxideSSE.git@v0.1.15
```

### Option 2. Clone and install in editable mode

Use this option if you want to modify the source code, add example files, or develop new functions.

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
pip install -e .
```

Check the installation:

```bash
python -c "import oxidesse as sse; print(sse.__version__)"
```

### Option 3. Reproducible conda environment

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
conda env create -f environment.yml
conda activate oxidesse
pip install -e .
```

---

## 🧠 SevenNet requirement

OxideSSE uses machine-learning interatomic potentials for **geometry optimization** and **stability analysis**. The default workflow is based on **SevenNet**, but SevenNet is **not installed automatically** by OxideSSE because GPU, CUDA, PyTorch, and model compatibility can depend strongly on the user's computing environment.

Install SevenNet separately before using SevenNet-based calculators:

- https://github.com/MDIL-SNU/SevenNet
- https://sevennet.readthedocs.io/en/latest/

### Supported calculator names

| Calculator name | Notes |
|---|---|
| `7net-0` | Default calculator. |
| `7net-l3i5` | SevenNet model name passed to `SevenNetCalculator`. |
| `7net-omat` | SevenNet model name passed to `SevenNetCalculator`. |
| `7net-mf-ompa` | Modal model. OxideSSE currently uses `modal="mpa"` by default; edit the calculator settings directly if another modal is needed. |
| `7net-omni` | Modal model. A modal must be specified by the user; edit the calculator settings or pass a manually constructed ASE calculator. |

Default:

```python
calculator="7net-0"
```

Advanced users can also pass any already-created ASE calculator object directly to OxideSSE functions.

---

## 🧪 Tested environment

OxideSSE was tested in the following environment:

| Component | Tested version |
|---|---:|
| Python | 3.12.13 |
| ASE | 3.28.0 |
| pymatgen | 2026.5.4 |
| mp-api | 0.46.1 |
| NumPy | 2.4.4 |
| pandas | 3.0.3 |
| SciPy | 1.17.1 |
| matplotlib | 3.10.9 |
| spglib | 2.7.0 |
| tqdm | 4.67.3 |
| fire | 0.7.1 |
| monty | 2026.2.18 |
| emmet-core | 0.86.4 |
| pydantic | 2.13.4 |
| pyarrow | 24.0.0 |
| SevenNet | 0.12.1, installed separately |

---

## 🔑 Materials Project API key

Formation energy and energy-above-hull calculations require access to the Materials Project API. Do not hard-code your API key in public scripts or notebooks.

Set it as an environment variable:

```bash
export MP_API_KEY="YOUR_MP_API_KEY"
```

Example for formation energy:

```python
form = sse.compute_binary_oxide_formation_energy(
    structure="Li7La3Zr2O12.cif",
)
```

Example for energy above hull:

```python
hull = sse.compute_energy_above_hull_mlp(
    structure="Li7La3Zr2O12.cif",
)
```

You can also pass the API key directly:

```python
form = sse.compute_binary_oxide_formation_energy(
    structure="Li7La3Zr2O12.cif",
    api_key="YOUR_MP_API_KEY",
)

hull = sse.compute_energy_above_hull_mlp(
    structure="Li7La3Zr2O12.cif",
    api_key="YOUR_MP_API_KEY",
)
```

---

# 🚀 Python API

## 1. ⚙️ `sse.optimize_structure()`

Optimize a CIF or POSCAR/CONTCAR-style structure file using an ASE-compatible calculator.

### Minimal example

```python
result = sse.optimize_structure(
    structure="input.cif",
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF, POSCAR, and CONTCAR-style files are supported. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `output_dir` | Directory where the optimized CIF and optimization log are saved. |
| `output_filename` | Name of the optimized structure file. If `None`, OxideSSE uses `opt_<input>.cif`. |
| `fmax` | Force convergence criterion in eV/Å. |
| `max_steps` | Maximum number of optimizer steps. |
| `cell_relax` | If `True`, relaxes both cell and atomic positions. If `False`, relaxes atomic positions only. |
| `device` | Device option passed to the calculator builder, for example `"auto"`, `"cpu"`, or `"cuda"`. |

### Output files

- Optimized CIF file, for example `optimized/opt_input.cif`
- Optimization log file, for example `optimized/input.opt.log`

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

Compute the formation energy of an oxide structure using binary oxide references. OxideSSE automatically balances the target composition against binary oxide references, fetches representative structures from Materials Project, recalculates their energies using the selected MLP calculator, and computes the formation energy.

Example composition balance:

```text
Li7La3Zr2O12 = 3.5 Li2O + 1.5 La2O3 + 2 ZrO2
```

### Minimal example

```python
result = sse.compute_binary_oxide_formation_energy(
    structure="input.cif",
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input oxide structure. CIF and POSCAR/CONTCAR-style files are supported. The structure must contain oxygen. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY` from the environment. |
| `output_csv` | CSV path where the result is saved. Set `None` to skip CSV output. |
| `custom_reference_formulas` | Optional mapping for unsupported elements, for example `{"Mg": "MgO"}`. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Output CSV

The CSV contains the target formula, target energy, formation energy, formation energy per atom, and the binary oxide references used in the calculation.

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

## 3. 📉 `sse.compute_energy_above_hull_mlp()`

Compute energy above hull using Materials Project entries as phase diagram references, while replacing MP energies with MLP-recalculated energies. OxideSSE uses MP structures and compositions, but the phase diagram is built from MLP energies.

### Minimal example

```python
result = sse.compute_energy_above_hull_mlp(
    structure="input.cif",
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF and POSCAR/CONTCAR-style files are supported. |
| `calculator` | ASE calculator object or supported calculator name. Default: `"7net-0"`. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY`. |
| `thermo_type` | MP thermo entry set. Supported values include `"R2SCAN"`, `"GGA_GGA+U"`, `"GGA_GGA+U_R2SCAN"`, `"GGAU"`, and `"GGAU_R2SCAN"`. |
| `output_csv` | CSV file for the hull result. |
| `cache_dir` | Directory where MLP-recalculated MP entry energies are cached. |
| `use_cache` | If `True`, reuse cached MLP energies for previously evaluated MP entries. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Output CSV

The CSV includes the formula, chemical system, thermo type, target MLP energy, energy above hull, number of MP entries used, failed entries, and cache path.

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

Analyze one LAMMPS simulation folder, convert the LAMMPS trajectory into structures, compute MSD curves, estimate diffusivity, perform error analysis, and save plots and CSV files.

### Minimal example

```python
result = sse.compute_diffusivity_from_lammps(
    simulation_dir="lammps_run",
    dump_file="dump.traj",
    data_file="structure.data",
    timestep_fs=1.0,
    step_skip=500,
    temperature=1000,
)
```

### Important variables

| Variable | Description |
|---|---|
| `simulation_dir` | Directory containing the LAMMPS simulation files. |
| `dump_file` | LAMMPS dump trajectory file. Absolute paths and paths relative to `simulation_dir` are supported. |
| `data_file` | LAMMPS data file. Absolute paths and paths relative to `simulation_dir` are supported. |
| `timestep_fs` | LAMMPS MD timestep in femtoseconds. For example, `1.0` means one LAMMPS step is 1 fs. |
| `step_skip` | Dump interval in MD steps. If frames are saved every 500 steps and `timestep_fs=1.0`, the frame spacing is 500 fs = 0.5 ps. |
| `temperature` | MD temperature in K. |
| `species` | Primary mobile species. Default: `"Li"`. |
| `output_dir` | Directory for MSD plot and CSV files. |
| `save_poscars` | If `True`, save converted trajectory frames as POSCAR files. |
| `save_plot` | If `True`, save the MSD plot. |
| `output_csv` | Summary CSV path. |
| `msd_csv` | Time-resolved MSD CSV path. |
| `plot_other_species` | If `True`, include non-primary species in the MSD plot. |
| `oxidized_species` | Optional oxidation-state mapping for conductivity-related calculations, for example `{"Li": "Li+"}`. |

### Output files

- `msd_<species>_<temperature>K.png`
- `diffusivity_summary.csv`
- `msd_by_species.csv`
- Optional converted POSCAR files if `save_poscars=True`

### Result attributes

```python
result.temperature
result.timestep_fs
result.step_skip
result.primary_species
result.plot_path
result.summary_csv_path
result.msd_csv_path
result.species_results
```

---

## 5. 🔥 `sse.plot_arrhenius()`

Fit Arrhenius behavior from temperature-dependent diffusivity data, save an Arrhenius plot, extrapolate diffusivity to 298 K, and optionally estimate 298 K ionic conductivity.

### Minimal example

```python
data = {
    600: {"diffusivity": 1.2e-7, "std": 0.2e-7},
    700: {"diffusivity": 5.0e-7, "std": 0.6e-7},
    800: {"diffusivity": 1.8e-6, "std": 0.2e-6},
}

result = sse.plot_arrhenius(data=data)
```

### Important variables

| Variable | Description |
|---|---|
| `data` | Dictionary containing temperature, diffusivity, and standard deviation. Required format: `{T: {"diffusivity": D, "std": std}}`. |
| `output_dir` | Directory where the plot and CSV are saved. |
| `output_filename` | Arrhenius plot filename. |
| `output_csv` | Summary CSV filename. Set `None` to skip CSV output. |
| `name` | Label for the measured diffusivity data. |
| `extrapolate_298K` | If `True`, plot the fitted line to 298 K and report 298 K diffusivity. |
| `structure` | Structure file used to compute 298 K ionic conductivity via the Nernst-Einstein relation. |
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

# 🖥️ Command-line interface

OxideSSE provides a Fire-based command-line interface after installation:

```bash
oxidesse --help
```

## Geometry optimization

```bash
oxidesse optimize-structure \
  --structure input.cif
```

Common options:

- `--structure`: input CIF/POSCAR path
- `--calculator`: calculator name, default `7net-0`
- `--output_dir`: output directory
- `--output_filename`: optimized CIF filename
- `--fmax`: force convergence criterion
- `--max_steps`: maximum optimization steps
- `--cell_relax`: whether to relax the cell
- `--device`: calculator device option

Batch optimize structures in a directory:

```bash
oxidesse optimize-structures \
  --input_dir cifs
```

## Formation energy

```bash
oxidesse formation-energy \
  --structure input.cif
```

Common options:

- `--structure`: input oxide structure
- `--calculator`: calculator name, default `7net-0`
- `--api_key`: Materials Project API key, optional if `MP_API_KEY` is set
- `--output_csv`: CSV output path
- `--device`: calculator device option
- `--append`: append to existing CSV

## Energy above hull

```bash
oxidesse energy-above-hull \
  --structure input.cif
```

Common options:

- `--structure`: input structure
- `--calculator`: calculator name, default `7net-0`
- `--api_key`: Materials Project API key, optional if `MP_API_KEY` is set
- `--thermo_type`: `R2SCAN`, `GGA_GGA+U`, `GGA_GGA+U_R2SCAN`, `GGAU`, or `GGAU_R2SCAN`
- `--output_csv`: CSV output path
- `--cache_dir`: cache directory
- `--use_cache`: reuse cached MLP entry energies
- `--append`: append to existing CSV

## Diffusivity from LAMMPS

```bash
oxidesse diffusivity \
  --simulation_dir lammps_run \
  --dump_file dump.traj \
  --data_file structure.data \
  --timestep_fs 1.0 \
  --step_skip 500 \
  --temperature 1000
```

Common options:

- `--simulation_dir`: LAMMPS simulation directory
- `--dump_file`: LAMMPS dump file
- `--data_file`: LAMMPS data file
- `--timestep_fs`: MD timestep in fs
- `--step_skip`: dump interval in MD steps
- `--temperature`: temperature in K
- `--species`: primary mobile species, default `Li`
- `--output_dir`: output directory
- `--save_poscars`: save converted POSCAR frames
- `--save_plot`: save MSD plot
- `--plot_other_species`: include non-primary species in plot

## Arrhenius fitting

For command-line use, pass data as a Python-style dictionary string:

```bash
oxidesse arrhenius \
  --data '{600: {"diffusivity": 1.2e-7, "std": 0.2e-7}, 700: {"diffusivity": 5.0e-7, "std": 0.6e-7}, 800: {"diffusivity": 1.8e-6, "std": 0.2e-6}}'
```

Common options:

- `--data`: temperature-diffusivity dictionary
- `--name`: dataset label
- `--output_dir`: output directory
- `--output_filename`: plot filename
- `--output_csv`: CSV filename
- `--extrapolate_298K`: whether to extrapolate to 298 K
- `--structure`: structure used for conductivity conversion
- `--specie`: oxidized mobile species, for example `Li+`
- `--conductivity_factor_298K`: optional precomputed conductivity conversion factor

---

# 📁 Example data

The `examples/` directory is intended for small files that users can run immediately after cloning the repository. Recommended contents are:

```text
examples/
├── structures/
│   ├── LLZO.cif
│   └── POSCAR_LLZO
└── lammps_run/
    ├── dump.traj
    ├── structure.data
    ├── log.lammps
    └── README.md
```

Recommended files to add:

- A small CIF file for geometry optimization and stability examples.
- A small POSCAR file for users who prefer VASP-style structure input.
- One short LAMMPS dump trajectory, ideally only a small number of frames.
- The matching LAMMPS data file for the dump trajectory.
- A short `README.md` inside `examples/lammps_run/` describing `timestep_fs`, `step_skip`, temperature, unit style, and species mapping.
- Optional small CSV file with temperature, diffusivity, and standard deviation values for Arrhenius testing.

Avoid uploading large production trajectories, unpublished screening datasets, cache directories, or raw Materials Project data.

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
├── tests/
├── environment.yml
├── pyproject.toml
└── README.md
```

---

# 📚 Citation

If you use OxideSSE in your research, please cite:

```bibtex
@misc{jeon2026unlocking,
  title = {Unlocking the Li7M3X2O12 Garnet Electrolyte Landscape with Universal Machine Learning Interatomic Potentials},
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
