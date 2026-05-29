import oxidesse as sse


def test_arrhenius_fit_runs():
    data = {
        600: {"diffusivity": 1e-7, "std": 1e-8},
        700: {"diffusivity": 5e-7, "std": 1e-8},
        800: {"diffusivity": 1e-6, "std": 1e-8},
    }
    result = sse.fit_arrhenius(data)
    assert result.activation_energy_eV > 0
    assert result.diffusivity_298K_cm2_s is not None
