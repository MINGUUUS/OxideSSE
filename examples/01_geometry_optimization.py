import oxidesse as sse

result = sse.optimize_structure(
    structure="input.cif",
    calculator="7net-0",
    output_dir="optimized",
    fmax=0.01,
    cell_relax=True,
)

print(result.converged)
print(result.energy_per_atom)
print(result.output_path)
