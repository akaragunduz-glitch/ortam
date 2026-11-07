"""ORTAM model contract.

Every module (groundwater, vadose zone, air, risk, ...) implements this
contract, so that every module can:
  1. run standalone (from the command line, from Python, or from the UI), and
  2. be used as one link in an integrated workflow (e.g. the risk chain).

Fortran analogy: this class plays the role of a shared INTERFACE block.
Every module provides the same three procedures: validate (input check),
run (computation) and summary (human-readable result).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import ClassVar

from pydantic import BaseModel


class ModelInput(BaseModel):
    """Base class for all module inputs. Every field documents its unit."""

    model_config = {"extra": "forbid"}  # An undeclared parameter is an error


class ModelOutput(BaseModel):
    """Base class for all module outputs."""

    warnings: list[str] = []


class Module(ABC):
    """An ORTAM module that can run standalone."""

    name: ClassVar[str]             # CLI name, e.g. "groundwater.three-point"
    title: ClassVar[str]            # Human-readable name
    input_type: ClassVar[type[ModelInput]]
    output_type: ClassVar[type[ModelOutput]]

    def validate(self, data: dict) -> ModelInput:
        """Validate raw data (YAML, JSON, form) and convert it to typed input."""
        return self.input_type.model_validate(data)

    @abstractmethod
    def run(self, inputs: ModelInput) -> ModelOutput:
        """Run the model. Must not depend on the database or other modules."""

    def summary(self, output: ModelOutput) -> str:
        """Short human-readable summary of the result (for the CLI)."""
        return output.model_dump_json(indent=2)
