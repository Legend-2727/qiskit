# CAES: Calibration-Aware Error-weighted SABRE

## Overview

CAES (Calibration-Aware Error-weighted SABRE) extends Qiskit's SABRE routing algorithm with error-awareness when backend calibration data is available.

## Usage

```python
from qiskit.transpiler.passes import SabreSwap

# Basic error-aware routing
sabre = SabreSwap(target, heuristic='error_aware')

# Fine-tuned parameters
sabre = SabreSwap(target, heuristic='error_aware', lambda_=0.6, omega=1.0)
```

## Parameters

- `lambda_` (0.0-1.0): Balance between lookahead and error minimization
- `omega` (float): Readout error penalty scaling  
- `distance_policy`: 'auto', 'hop', or 'error_weighted'

## How It Works

CAES converts error probabilities to distance weights using `w = -3 * ln(1 - p)` and integrates them into SABRE's routing decisions. When no calibration data is available, it falls back to standard SABRE behavior.

## Backward Compatibility

All existing SabreSwap usage continues unchanged. The feature is opt-in only through `heuristic='error_aware'`.
