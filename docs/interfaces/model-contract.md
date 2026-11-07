# Module contract

Every ORTAM module subclasses `Module` in `ortam_core/model.py`.

| Member | Description |
|--------|-------------|
| `name` | CLI name, `<module>.<tool>` (e.g. `groundwater.three-point`, `groundwater.modflow`) |
| `title` | Human-readable name |
| `input_type` | Pydantic class; every field documents its unit; undeclared fields are errors |
| `output_type` | Pydantic class; every output has a `warnings` field |
| `validate(data)` | Converts raw data to typed input; invalid input is rejected here |
| `run(inputs)` | Performs the computation; no dependency on the database or other modules |
| `summary(output)` | Short human-readable summary for the CLI |

## Checklist for a new module

- [ ] `modules/<module>/<tool>.py` — class, input and output types
- [ ] `modules/<module>/examples/<tool>.yaml` — a working example input
- [ ] `modules/<module>/tests/test_<tool>.py` — reference tests
- [ ] `modules/registry.py` — one-line registration
- [ ] `docs/methods/<module>.md` — equations, parameters, references

Reference implementation: `modules/groundwater/three_point.py`.
