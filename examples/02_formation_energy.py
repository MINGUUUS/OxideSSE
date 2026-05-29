import oxidesse as sse

result = sse.compute_binary_oxide_formation_energy(
    structure="input.cif",
    calculator="7net-0",
    output_csv="formation_energy.csv",
)

print(result.formula)
print(result.formation_energy_per_atom)
for ref in result.references:
    print(ref.formula, ref.coefficient, ref.energy_per_formula)
