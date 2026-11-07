"""Defaults for the air emission-factor tools."""

from __future__ import annotations

from functools import lru_cache

import yaml

from ortam_core.chemicals.rsl_import import DATA


@lru_cache(maxsize=None)
def defaults(version: str = "2024-11") -> dict[str, dict[str, float]]:
    f = yaml.safe_load((DATA / version / "air_emission_factors.yaml").read_text(encoding="utf-8"))
    return {g: {k: float(v["value"]) for k, v in fields.items()} for g, fields in f.items() if isinstance(fields, dict)}
