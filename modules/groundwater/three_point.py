"""Three-point method: hydraulic gradient and groundwater flow direction.

Method: solve for the plane h(x, y) = a*x + b*y + c through the hydraulic
heads of three wells. Hydraulic gradient i = sqrt(a^2 + b^2); flow is in the
direction of -grad(h).
Source: classic three-point problem (e.g. Fetter, Applied Hydrogeology); the
plane equation is analytic geometry. Verification: hand calculation (see tests).

This tool does not require MODFLOW; it gives a quick first assessment from
field data. The numerical flow model (MODFLOW 6) will be added to the same
module later as `groundwater.modflow`.
"""

from __future__ import annotations

import math

import numpy as np
from pydantic import Field

from ortam_core.model import ModelInput, ModelOutput, Module


class Well(ModelInput):
    name: str
    x_m: float = Field(description="Easting, m (metric projection, e.g. UTM)")
    y_m: float = Field(description="Northing, m")
    head_m: float = Field(description="Hydraulic head, m (top of casing elevation - depth to water)")


class ThreePointInput(ModelInput):
    wells: list[Well] = Field(min_length=3, max_length=3)
    K_m_d: float | None = Field(None, gt=0, description="Hydraulic conductivity, m/d")
    n_e: float | None = Field(None, gt=0, le=1, description="Effective porosity, -")


class ThreePointOutput(ModelOutput):
    gradient: float = Field(description="Hydraulic gradient i, m/m")
    flow_azimuth_deg: float = Field(description="Flow direction, degrees clockwise from north")
    darcy_flux_m_d: float | None = None
    seepage_velocity_m_d: float | None = None


class ThreePointModule(Module):
    name = "groundwater.three-point"
    title = "Groundwater: three-point flow direction"
    input_type = ThreePointInput
    output_type = ThreePointOutput

    def run(self, inputs: ThreePointInput) -> ThreePointOutput:
        # Fortran equivalent: solve the 3x3 linear system A * [a, b, c] = h
        A = np.array([[w.x_m, w.y_m, 1.0] for w in inputs.wells])
        h = np.array([w.head_m for w in inputs.wells])
        if abs(np.linalg.det(A)) < 1e-9:
            raise ValueError("Wells are collinear; the three-point solution is undefined.")
        a, b, _ = np.linalg.solve(A, h)

        gradient = math.hypot(a, b)
        # Flow along -grad(h): east component -a, north component -b
        azimuth = math.degrees(math.atan2(-a, -b)) % 360.0

        warnings = []
        if gradient < 1e-4:
            warnings.append("Very low gradient (< 1e-4); flow direction may be uncertain.")

        q = v = None
        if inputs.K_m_d is not None:
            q = inputs.K_m_d * gradient          # Darcy's law, q = K * i
            if inputs.n_e is not None:
                v = q / inputs.n_e               # Average linear (seepage) velocity

        return ThreePointOutput(gradient=gradient, flow_azimuth_deg=azimuth,
                                darcy_flux_m_d=q, seepage_velocity_m_d=v,
                                warnings=warnings)

    def summary(self, o: ThreePointOutput) -> str:
        lines = [
            f"Hydraulic gradient   : {o.gradient:.5f} m/m",
            f"Flow direction       : {o.flow_azimuth_deg:.1f} deg (clockwise from north)",
        ]
        if o.darcy_flux_m_d is not None:
            lines.append(f"Darcy flux           : {o.darcy_flux_m_d:.4g} m/d")
        if o.seepage_velocity_m_d is not None:
            lines.append(f"Seepage velocity     : {o.seepage_velocity_m_d:.4g} m/d")
        lines += [f"WARNING: {w}" for w in o.warnings]
        return "\n".join(lines)
