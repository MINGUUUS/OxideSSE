# LAMMPS trajectory example

Add a small LAMMPS test run here for diffusivity analysis.

Suggested files:

- `dump.traj`: short LAMMPS dump trajectory
- `structure.data`: corresponding LAMMPS data file
- `log.lammps`: optional simulation log
- `README.md`: timestep, dump interval, temperature, and units used in the run

The folder can be used with:

```python
import oxidesse as sse
sse.compute_diffusivity_from_lammps(
    simulation_dir="examples/lammps_run",
    dump_file="dump.traj",
    data_file="structure.data",
    timestep_fs=1.0,
    step_skip=500,
    temperature=1000,
)
```
