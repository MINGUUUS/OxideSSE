import oxidesse as sse

result = sse.compute_diffusivity_from_lammps(
    simulation_dir="lammps_run",
    dump_file="dump.traj",
    data_file="structure.data",
    timestep_fs=1.0,
    step_skip=500,
    temperature=1000,
    species="Li",
    output_dir="diffusivity_results",
    oxidized_species={"Li": "Li+"},
)

print(result.diffusivity)
print(result.diffusivity_std)
print(result.summary_csv_path)
