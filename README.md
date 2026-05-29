# 🧪 OxideSSE

**OxideSSE** is a Python toolkit for oxide solid-state electrolyte research, with a focus on LLZO-derived garnet-type oxide materials. It provides reusable workflows for structure relaxation, stability screening, LAMMPS trajectory analysis, diffusivity calculation, and Arrhenius extrapolation.

The package is designed to be used as:

```python
import oxidesse as sse
```

## ✨ Main features

- ⚙️ Geometry optimization using ASE-compatible MLP calculators.
- 🧱 Binary-oxide-referenced formation energy with MLP-recalculated reference energies.
- 📉 Energy above hull using Materials Project entries with MLP-recalculated entry energies.
- 🚶 Li-ion MSD and diffusivity analysis from a single LAMMPS simulation folder.
- 🔥 Arrhenius fitting, 298 K diffusivity extrapolation, and ionic conductivity estimation.

## 📦 Installation

OxideSSE can be installed directly from GitHub or cloned for development.

### Option 1. Install directly from GitHub

Use this option if you only want to use the package.

```bash
pip install git+https://github.com/MINGUUUS/OxideSSE.git
```

To install a specific release tag:

```bash
pip install git+https://github.com/MINGUUUS/OxideSSE.git@v0.1.14
```

### Option 2. Clone and install in editable mode

Use this option if you want to modify the source code, run examples, or develop new functions.

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
pip install -e .
```

After installation, check that OxideSSE can be imported:

```bash
python -c "import oxidesse as sse; print(sse.__version__)"
```

### Reproducible conda environment

For better reproducibility, you can create a conda environment from the provided `environment.yml` file:

```bash
git clone https://github.com/MINGUUUS/OxideSSE.git
cd OxideSSE
conda env create -f environment.yml
conda activate oxidesse
pip install -e .
```

The `environment.yml` file records the versions used during package testing. If SevenNet installation fails because of CUDA/PyTorch compatibility, install SevenNet manually by following the official documentation.

### SevenNet requirement

This project uses **SevenNet** as the default MLP backend. SevenNet should be installed in the same Python environment before using:

```python
calculator="7net-0"
calculator="7net-d3"
calculator="7net-mf-ompa"
```

Please follow the official SevenNet installation instructions:

- https://github.com/MDIL-SNU/SevenNet
- https://sevennet.readthedocs.io/en/latest/

## 🧪 Tested environment

OxideSSE v0.1.14 was tested in the following environment:

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
| SevenNet | 0.12.1 |

`pyproject.toml` uses compatible version ranges for normal installation, while `environment.yml` provides a more reproducible tested environment.

## 🔑 Materials Project API key

Formation energy and energy-above-hull calculations require access to the Materials Project API. Do not hard-code your API key in public scripts or notebooks.

Set it as an environment variable:

```bash
export MP_API_KEY="YOUR_MP_API_KEY"
```

Then call the functions normally:

```python
form = sse.compute_binary_oxide_formation_energy(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
)

hull = sse.compute_energy_above_hull_mlp(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
)
```

You can also pass the API key directly:

```python
form = sse.compute_binary_oxide_formation_energy(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
    api_key="YOUR_MP_API_KEY",
)

hull = sse.compute_energy_above_hull_mlp(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
    api_key="YOUR_MP_API_KEY",
)
```

## 🚀 Quick example

```python
import oxidesse as sse

opt = sse.optimize_structure(
    structure="input.cif",
    calculator="7net-0",
    output_dir="optimized",
)

form = sse.compute_binary_oxide_formation_energy(
    structure=opt.output_path,
    calculator="7net-0",
    output_csv="formation_energy.csv",
)

hull = sse.compute_energy_above_hull_mlp(
    structure=opt.output_path,
    calculator="7net-0",
    thermo_type="R2SCAN",
    output_csv="hull_energy.csv",
)
```

---

# 🧩 Python API

## 1. ⚙️ `sse.optimize_structure()`

Optimize a CIF or POSCAR/CONTCAR-style structure file using an ASE-compatible calculator.

### Input

```python
result = sse.optimize_structure(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
    output_dir="optimized",
    output_filename=None,
    fmax=0.01,
    max_steps=500,
    cell_relax=True,
    device="auto",
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF, POSCAR, and CONTCAR-style files are supported. |
| `calculator` | ASE calculator object or supported calculator name such as `"7net-0"`, `"7net-d3"`, or `"7net-mf-ompa"`. |
| `output_dir` | Directory where the optimized CIF and optimization log are saved. |
| `output_filename` | Name of the optimized structure file. If `None`, OxideSSE uses `opt_<input>.cif`. |
| `fmax` | Force convergence criterion in eV/Å. |
| `max_steps` | Maximum number of optimizer steps. |
| `cell_relax` | If `True`, relaxes both cell and atomic positions. If `False`, relaxes atomic positions only. |
| `device` | Device option passed to the calculator builder, for example `"auto"`, `"cpu"`, or `"cuda"`. |

### Output files

- Optimized CIF file, for example `optimized/opt_Li7La3Zr2O12.cif`
- Optimization log file, for example `optimized/Li7La3Zr2O12.opt.log`

### Result attributes

```python
print(result.converged)       # True/False
print(result.energy)          # Final total energy in eV
print(result.energy_per_atom) # Final energy per atom in eV/atom
print(result.max_force)       # Final maximum force in eV/Å
print(result.output_path)     # Path to optimized CIF
print(result.log_path)        # Path to log file
print(result.message)         # Status message
```

---

## 2. 🧱 `sse.compute_binary_oxide_formation_energy()`

Compute the formation energy of an oxide structure using binary oxide references. OxideSSE automatically balances the target composition against binary oxide references, fetches representative structures from Materials Project, recalculates their energies using the selected MLP calculator, and computes the formation energy.

For example:

```text
Li7La3Zr2O12 = 3.5 Li2O + 1.5 La2O3 + 2 ZrO2
```

No reference energy is hard-coded. Even `Li2O` is fetched and recalculated using the MLP calculator.

### Input

```python
result = sse.compute_binary_oxide_formation_energy(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
    api_key=None,
    output_csv="formation_energy.csv",
    custom_reference_formulas=None,
    device="auto",
    append=True,
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input oxide structure. CIF and POSCAR/CONTCAR-style files are supported. The structure must contain oxygen. |
| `calculator` | ASE calculator object or supported calculator name. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY` from the environment. |
| `output_csv` | CSV path where the result is saved. Set `None` to skip CSV output. |
| `custom_reference_formulas` | Optional mapping for unsupported elements, for example `{"Mg": "MgO"}`. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Reference rules

- `Li` → `Li2O`
- Elements in `B3_ELEMENTS` → `X2O3`
- Elements in `C4_ELEMENTS` → `XO2`

If an element is not covered by the default rules, use `custom_reference_formulas`.

### Output CSV

The CSV contains the target formula, target energy, formation energy, formation energy per atom, and the binary oxide references used in the calculation.

### Result attributes

```python
print(result.formula)
print(result.target_energy)                  # Target MLP total energy in eV
print(result.target_energy_per_atom)         # Target MLP energy per atom
print(result.formation_energy)               # Formation energy in eV per input cell/formula
print(result.formation_energy_per_atom)      # Formation energy in eV/atom
print(result.csv_path)                       # Output CSV path

for ref in result.references:
    print(ref.element)
    print(ref.formula)
    print(ref.coefficient)
    print(ref.material_id)
    print(ref.energy_per_formula)
```

---

## 3. 📉 `sse.compute_energy_above_hull_mlp()`

Compute energy above hull using Materials Project entries as phase diagram references, while replacing MP energies with MLP-recalculated energies.

This means OxideSSE uses MP structures and compositions, but the phase diagram is built from MLP energies.

### Input

```python
result = sse.compute_energy_above_hull_mlp(
    structure="Li7La3Zr2O12.cif",
    calculator="7net-0",
    api_key=None,
    thermo_type="R2SCAN",
    output_csv="hull_energy.csv",
    cache_dir=".oxidesse_cache",
    use_cache=True,
    device="auto",
    append=True,
)
```

### Important variables

| Variable | Description |
|---|---|
| `structure` | Input structure file. CIF and POSCAR/CONTCAR-style files are supported. |
| `calculator` | ASE calculator object or supported calculator name. |
| `api_key` | Materials Project API key. If `None`, OxideSSE reads `MP_API_KEY`. |
| `thermo_type` | MP thermo entry set. Supported values include `"R2SCAN"`, `"GGA_GGA+U"`, `"GGA_GGA+U_R2SCAN"`, `"GGAU"`, and `"GGAU_R2SCAN"`. |
| `output_csv` | CSV file for the hull result. |
| `cache_dir` | Directory where MLP-recalculated MP entry energies are cached. |
| `use_cache` | If `True`, reuse cached MLP energies for previously evaluated MP entries. |
| `device` | Calculator device option. |
| `append` | If `True`, append to an existing CSV file. |

### Output CSV

The CSV includes:

- formula
- chemical system
- thermo type
- target MLP energy
- target MLP energy per atom
- energy above hull
- number of MP entries used
- number of failed entries
- cache path

### Result attributes

```python
print(result.formula)
print(result.chemical_system)
print(result.thermo_type)
print(result.target_energy)
print(result.target_energy_per_atom)
print(result.energy_above_hull)
print(result.n_entries)
print(result.failed_entries)
print(result.cache_path)
print(result.csv_path)
```

---

## 4. 🚶 `sse.compute_diffusivity_from_lammps()`

Analyze one LAMMPS simulation folder, convert the LAMMPS trajectory into structures, compute MSD curves, estimate diffusivity, perform error analysis, and save plots and CSV files.

### Input

```python
result = sse.compute_diffusivity_from_lammps(
    simulation_dir="lammps_run",
    dump_file="dump.traj",
    data_file="structure.data",
    timestep_fs=1.0,
    step_skip=500,
    temperature=1000,
    species="Li",
    output_dir="diffusivity_results",
    save_poscars=False,
    save_plot=True,
    output_csv="diffusivity_summary.csv",
    msd_csv="msd_by_species.csv",
    plot_other_species=True,
    oxidized_species={"Li": "Li+"},
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
| `species` | Primary mobile species, usually `"Li"`. |
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

The MSD plot is always shown in picoseconds. For Li, the plot includes the linear fitting region used for the diffusivity calculation, and the title includes the fitted diffusivity in cm²/s.

### Result attributes

```python
print(result.temperature)
print(result.timestep_fs)
print(result.step_skip)
print(result.primary_species)
print(result.plot_path)
print(result.summary_csv_path)
print(result.msd_csv_path)

for item in result.species_results:
    print(item.species)
    print(item.diffusivity)
    print(item.diffusivity_std)
    print(item.diffusivity_relative_std)
    print(item.diffusivity_components)
    print(item.n_jump)
    print(item.n_jump_component)
    print(item.fit_start_time_ps)
    print(item.fit_end_time_ps)
```

---

## 5. 🔥 `sse.plot_arrhenius()`

Fit Arrhenius behavior from temperature-dependent diffusivity data, save an Arrhenius plot, extrapolate diffusivity to 298 K, and optionally estimate 298 K ionic conductivity.

### Input

```python
data = {
    600: {"diffusivity": 1.2e-7, "std": 0.2e-7},
    700: {"diffusivity": 5.0e-7, "std": 0.6e-7},
    800: {"diffusivity": 1.8e-6, "std": 0.2e-6},
}

result = sse.plot_arrhenius(
    data=data,
    output_dir="arrhenius_results",
    output_filename="arrhenius_plot.png",
    output_csv="arrhenius_summary.csv",
    name="LLZO",
    extrapolate_298K=True,
    structure="Li7La3Zr2O12.cif",
    specie="Li+",
    conductivity_factor_298K=None,
)
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

The 298 K extrapolated point is plotted with an error range when uncertainty can be estimated from the Arrhenius fit. It is not given a separate legend label.

### Result attributes

```python
print(result.activation_energy_eV)
print(result.activation_energy_error_eV)
print(result.diffusivity_298K_cm2_s)
print(result.diffusivity_298K_std_cm2_s)
print(result.diffusivity_298K_range_cm2_s)
print(result.conductivity_298K_mS_cm)
print(result.conductivity_298K_range_mS_cm)
print(result.plot_path)
print(result.csv_path)
```

---

# 🖥️ Command-line interface

OxideSSE provides a Fire-based command-line interface after installation:

```bash
oxidesse --help
```

## Geometry optimization

Optimize one structure:

```bash
oxidesse optimize-structure \
  --structure input.cif \
  --calculator 7net-0 \
  --output_dir optimized \
  --fmax 0.01 \
  --max_steps 500 \
  --cell_relax True \
  --device auto
```

Useful variables:

- `--structure`: input CIF/POSCAR path
- `--calculator`: calculator name
- `--output_dir`: output directory
- `--output_filename`: optimized CIF filename
- `--fmax`: force convergence criterion
- `--max_steps`: maximum optimization steps
- `--cell_relax`: whether to relax the cell
- `--device`: calculator device option

Batch optimize structures in a directory:

```bash
oxidesse optimize-structures \
  --input_dir cifs \
  --pattern "*.cif" \
  --output_dir optimized \
  --calculator 7net-0
```

## Formation energy

```bash
oxidesse formation-energy \
  --structure input.cif \
  --calculator 7net-0 \
  --output_csv formation_energy.csv \
  --device auto
```

Useful variables:

- `--structure`: input oxide structure
- `--calculator`: calculator name
- `--api_key`: Materials Project API key, optional if `MP_API_KEY` is set
- `--output_csv`: CSV output path
- `--device`: calculator device option
- `--append`: append to existing CSV

## Energy above hull

```bash
oxidesse energy-above-hull \
  --structure input.cif \
  --calculator 7net-0 \
  --thermo_type R2SCAN \
  --output_csv hull_energy.csv \
  --cache_dir .oxidesse_cache \
  --use_cache True
```

Useful variables:

- `--structure`: input structure
- `--calculator`: calculator name
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
  --temperature 1000 \
  --species Li \
  --output_dir diffusivity_results
```

Useful variables:

- `--simulation_dir`: LAMMPS simulation directory
- `--dump_file`: LAMMPS dump file
- `--data_file`: LAMMPS data file
- `--timestep_fs`: MD timestep in fs
- `--step_skip`: dump interval in MD steps
- `--temperature`: temperature in K
- `--species`: primary mobile species
- `--output_dir`: output directory
- `--save_poscars`: save converted POSCAR frames
- `--save_plot`: save MSD plot
- `--plot_other_species`: include non-primary species in plot

## Arrhenius fitting

For command-line use, pass data as a Python-style dictionary string:

```bash
oxidesse arrhenius \
  --data '{600: {"diffusivity": 1.2e-7, "std": 0.2e-7}, 700: {"diffusivity": 5.0e-7, "std": 0.6e-7}, 800: {"diffusivity": 1.8e-6, "std": 0.2e-6}}' \
  --name LLZO \
  --output_dir arrhenius_results \
  --extrapolate_298K True \
  --structure Li7La3Zr2O12.cif \
  --specie Li+
```

Useful variables:

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
├── pyproject.toml
└── README.md
```
