import oxidesse as sse

data = {
    600: {"diffusivity": 1.2e-7, "std": 0.2e-7},
    700: {"diffusivity": 5.0e-7, "std": 0.6e-7},
    800: {"diffusivity": 1.8e-6, "std": 0.2e-6},
}

result = sse.plot_arrhenius(
    data,
    name="example",
    output_dir="arrhenius_results",
    extrapolate_298K=True,
    structure=None,  # pass CIF/POSCAR to also compute ionic conductivity
    specie="Li+",
)

print(result.activation_energy_eV)
print(result.diffusivity_298K_cm2_s)
print(result.conductivity_298K_mS_cm)
