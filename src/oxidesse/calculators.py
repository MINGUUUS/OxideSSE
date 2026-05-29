"""Calculator factory utilities."""
from __future__ import annotations

from typing import Any


def get_calculator(calculator: str | Any = "7net-0", device: str = "auto", **kwargs: Any):
    """Return an ASE calculator.

    Parameters
    ----------
    calculator
        Calculator name or an already-created ASE calculator. Supported names are
        ``"7net-0"``, ``"7net-d3"``, and ``"7net-mf-ompa"``.
    device
        Device string passed to SevenNet where applicable.
    **kwargs
        Additional keyword arguments forwarded to the calculator constructor.
    """
    if not isinstance(calculator, str):
        return calculator

    name = calculator.lower()
    try:
        from sevenn.calculator import SevenNetCalculator, SevenNetD3Calculator
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "SevenNet is required for calculator='7net-*'. Install SevenNet following https://github.com/MDIL-SNU/SevenNet or https://sevennet.readthedocs.io/en/latest/, or pass "
            "an already-created ASE calculator object."
        ) from exc

    if name == "7net-0":
        return SevenNetCalculator("7net-0", device=device, **kwargs)
    if name == "7net-d3":
        model = kwargs.pop("model", "7net-0")
        return SevenNetD3Calculator(model, device=device, **kwargs)
    if name == "7net-mf-ompa":
        return SevenNetCalculator("7net-mf-ompa", modal="mpa", **kwargs)

    raise ValueError(
        "Unsupported calculator. Use '7net-0', '7net-d3', '7net-mf-ompa', "
        "or pass an ASE calculator object."
    )


def calculator_label(calculator: str | Any) -> str:
    """Return a stable label for cache keys and CSV output."""
    if isinstance(calculator, str):
        return calculator
    return calculator.__class__.__name__
