# This code is part of Qiskit.
#
# (C) Copyright IBM 2025.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

"""Helper functions for extracting error maps from backends for CAES routing."""

import logging
from typing import Dict, Tuple, Optional, Set
from qiskit.transpiler.target import Target

LOG = logging.getLogger(__name__)


def extract_error_maps(target: Target) -> Tuple[Optional[Dict[Tuple[int, int], float]], Optional[Dict[int, float]]]:
    """
    Extract 2-qubit gate errors and readout errors from a Target.
    
    Args:
        target: Qiskit Target object containing backend properties
        
    Returns:
        Tuple of (edge_error_map, readout_error_map) where:
        - edge_error_map: Dict mapping (qubit_a, qubit_b) -> gate_error_rate
        - readout_error_map: Dict mapping qubit -> readout_error_rate  
        Returns (None, None) if no calibration data is available.
        
    Notes:
        - Gate errors are clamped to [0, 1) to avoid log(0) in distance calculations
        - Missing or invalid values are skipped
        - Returns None maps if no usable calibration data is found
    """
    edge_errors = {}
    readout_errors = {}
    
    try:
        # Extract 2-qubit gate errors by iterating through target operations
        # Look for common 2-qubit gates: cx, cnot, cz, iswap, etc.
        two_qubit_gates = {'cx', 'cnot', 'cz', 'iswap', 'ecr'}
        
        for gate_name in two_qubit_gates:
            if gate_name in target:
                # Access the gate directly from target
                gate_data = target[gate_name]
                for qargs, props in gate_data.items():
                    if len(qargs) == 2 and props is not None:
                        try:
                            if hasattr(props, 'error') and props.error is not None:
                                error_rate = float(props.error)
                                # Clamp to valid range - avoid p >= 1.0 which would cause log(0)
                                if 0.0 <= error_rate < 0.999:
                                    edge_key = tuple(sorted([qargs[0], qargs[1]]))  # Normalize edge order
                                    edge_errors[edge_key] = error_rate
                                else:
                                    LOG.debug(f"Skipping invalid error rate {error_rate} for {gate_name}{qargs}")
                        except (AttributeError, ValueError, TypeError) as e:
                            LOG.debug(f"Could not extract error for {gate_name}{qargs}: {e}")
                            continue

        # Extract readout errors from measure operations
        if 'measure' in target:
            measure_data = target['measure']
            for qargs, props in measure_data.items():
                if len(qargs) == 1 and props is not None:
                    try:
                        if hasattr(props, 'error') and props.error is not None:
                            readout_rate = float(props.error)
                            if 0.0 <= readout_rate < 1.0:
                                readout_errors[qargs[0]] = readout_rate
                    except (AttributeError, ValueError, TypeError) as e:
                        LOG.debug(f"Could not extract readout error for qubit {qargs[0]}: {e}")
                        continue
                        
    except Exception as e:
        LOG.warning(f"Error extracting calibration data from target: {e}")
        return None, None
    
    # Return results - None if no data found
    edge_map = edge_errors if edge_errors else None
    readout_map = readout_errors if readout_errors else None
    
    if edge_map or readout_map:
        LOG.debug(f"Extracted {len(edge_errors)} edge errors, {len(readout_errors)} readout errors")
        return edge_map, readout_map
    else:
        LOG.debug("No calibration data found in target")
        return None, None


def validate_error_maps(
    edge_errors: Optional[Dict[Tuple[int, int], float]], 
    readout_errors: Optional[Dict[int, float]],
    num_qubits: int
) -> Tuple[bool, str]:
    """
    Validate extracted error maps for consistency.
    
    Args:
        edge_errors: Edge error map to validate
        readout_errors: Readout error map to validate  
        num_qubits: Expected number of qubits in the system
        
    Returns:
        Tuple of (is_valid, error_message)
    """
    if edge_errors is None and readout_errors is None:
        return False, "No error data available"
        
    # Validate edge errors
    if edge_errors is not None:
        for (q1, q2), error_rate in edge_errors.items():
            if not (0 <= q1 < num_qubits and 0 <= q2 < num_qubits):
                return False, f"Edge ({q1}, {q2}) references invalid qubits for {num_qubits}-qubit system"
            if not (0.0 <= error_rate < 1.0):
                return False, f"Invalid error rate {error_rate} for edge ({q1}, {q2})"
                
    # Validate readout errors  
    if readout_errors is not None:
        for qubit, error_rate in readout_errors.items():
            if not (0 <= qubit < num_qubits):
                return False, f"Readout error for invalid qubit {qubit} in {num_qubits}-qubit system"
            if not (0.0 <= error_rate < 1.0):
                return False, f"Invalid readout error rate {error_rate} for qubit {qubit}"
                
    return True, ""


def should_use_error_weighting(
    heuristic: str,
    distance_policy: str, 
    edge_errors: Optional[Dict[Tuple[int, int], float]]
) -> bool:
    """
    Determine whether to use error-weighted distances based on policy and data availability.
    
    Args:
        heuristic: The routing heuristic being used
        distance_policy: User-specified distance policy ('auto', 'hop', 'error_weighted')
        edge_errors: Available edge error data (None if not available)
        
    Returns:
        True if error-weighted distances should be used, False for hop distances
    """
    if distance_policy == "hop":
        return False
    elif distance_policy == "error_weighted":
        return edge_errors is not None
    elif distance_policy == "auto":
        # Auto policy: use error weighting for error_aware heuristic if data is available
        return heuristic == "error_aware" and edge_errors is not None
    else:
        # Default to hop distances for unknown policies
        LOG.warning(f"Unknown distance_policy '{distance_policy}', defaulting to hop distances")
        return False