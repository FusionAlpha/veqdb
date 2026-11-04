"""Read and write VEQDB compact equilibrium records."""

from ._result import FitResult
from .compact import from_compact_json, to_compact_json

__all__ = ["FitResult", "from_compact_json", "to_compact_json"]
