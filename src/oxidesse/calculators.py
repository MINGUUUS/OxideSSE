"""Calculator factory utilities."""
from __future__ import annotations

from typing import Any


SUPPORTED_SEVENNET_CALCULATORS = (
    "7net-0",
    "7net-l3i5",
    "7net-omat",
    "7net-mf-ompa",
    "7net-omni",
    "7net-d3",
)


def get_calculator(calculator: str | Any = "7net-0", device: str = "auto", **kwargs: Any):
    """Return an ASE calculator.

    Parameters
    ----------
    calculator
        Calculator name or an already-created ASE calculator. Supported built-in
        SevenNet names are ``"7net-0"``, ``"7net-l3i5"``, ``"7net-omat"``,
        ``"7net-mf-ompa"``, ``"7net-omni"``, and ``"7net-d3"``.
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
            "SevenNet is required for calculator='7net-*'. Install SevenNet separately "
            "by following https://github.com/MDIL-SNU/SevenNet or "
            "https://sevennet.readthedocs.io/en/latest/, or pass an already-created "
            "ASE calculator object."
        ) from exc

    if name in {"7net-0", "7net-l3i5", "7net-omat"}:
        return SevenNetCalculator(name, device=device, **kwargs)

    if name == "7net-d3":
        model = kwargs.pop("model", "7net-0")
        return SevenNetD3Calculator(model, device=device, **kwargs)

    if name == "7net-mf-ompa":
        # 7net-mf-ompa is a modal model. The default modal is set to "mpa"
        # to match the original OxideSSE workflow, but advanced users may pass
        # another modal through keyword arguments when constructing calculators.
        modal = kwargs.pop("modal", "mpa")
        return SevenNetCalculator("7net-mf-ompa", modal=modal, device=device, **kwargs)

    if name == "7net-omni":
        # 7net-omni is modal-dependent. OxideSSE intentionally does not choose
        # a universal default modal because the correct modal depends on the
        # user's target domain. Pass modal="..." or create the calculator manually.
        modal = kwargs.pop("modal", None)
        if modal is None:
            raise ValueError(
                "calculator='7net-omni' requires a modal argument. Pass modal='...' "
                "through the Python API, or create a SevenNet calculator manually and "
                "pass it to OxideSSE."
            )
        return SevenNetCalculator("7net-omni", modal=modal, device=device, **kwargs)

    raise ValueError(
        "Unsupported calculator. Use one of "
        f"{SUPPORTED_SEVENNET_CALCULATORS}, or pass an ASE calculator object."
    )


def calculator_label(calculator: str | Any) -> str:
    """Return a stable label for cache keys and CSV output."""
    if isinstance(calculator, str):
        return calculator
    return calculator.__class__.__name__
