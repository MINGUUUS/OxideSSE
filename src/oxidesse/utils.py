"""General utilities for OxideSSE."""
from __future__ import annotations

from pathlib import Path


def make_directory(path: str | Path) -> Path:
    path = Path(path).expanduser().resolve()
    path.mkdir(parents=True, exist_ok=True)
    return path
