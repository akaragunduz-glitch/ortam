# ORTAM

**Open Risk, Transport And Modelling** — a modular environmental modelling
platform for soil, groundwater, soil gas, air and the food chain.

Every module (groundwater, vadose zone, vapor intrusion, air dispersion,
geochemistry, statistics, risk, remediation, ...) runs **standalone** or can
be connected to others through workflows.

```bash
uv run ortam list
uv run ortam run groundwater.three-point modules/groundwater/examples/three_point.yaml
```

- Scope: [docs/00-scope.md](docs/00-scope.md)
- Architecture: [docs/01-architecture.md](docs/01-architecture.md)
- Module contract: [docs/interfaces/model-contract.md](docs/interfaces/model-contract.md)
- Working rules for Claude: [CLAUDE.md](CLAUDE.md)

Status: Phase 0 — infrastructure setup.
