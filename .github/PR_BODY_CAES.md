## CAES: Calibration-Aware Error-Weighted SABRE (opt-in)

**What it does**: Adds an opt-in error-aware routing heuristic to SABRE:
- edge weights w = -3 ln(1-p) from Target 2q error rates
- optional readout penalty omega*(r_a + r_b)
- keeps SABRE lookahead via lambda

**Defaults unchanged**: If not enabled or if calibration is missing, routing is identical to baseline.

**How to enable**
```python
from qiskit.transpiler.passes.routing import SabreSwap
pm_pass = SabreSwap(target, heuristic="error_aware", lambda_=0.6, omega=1.0, seed=42)
```

**Tests/bench**: CAES unit tests pass; parity test passes; mini-bench CSV at repo root.

**Risks**: Dependent on calibration freshness; small runtime overhead due to APSP precompute (cached).
