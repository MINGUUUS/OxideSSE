import oxidesse as sse

result = sse.compute_energy_above_hull_mlp(
    structure="input.cif",
    calculator="7net-0",
    thermo_type="R2SCAN",
    output_csv="hull_energy.csv",
    cache_dir=".oxidesse_cache",
)

print(result.energy_above_hull)
print(result.n_entries)
