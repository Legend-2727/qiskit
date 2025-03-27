# Roadmap for Incorporating Error Rates into LightSABRE

This document outlines the modifications made to the Qiskit repository to incorporate error rates into the LightSABRE qubit mapping algorithm. These changes enhance circuit fidelity by dynamically adapting to hardware noise characteristics while maintaining scalability and runtime efficiency.

## Objective

The goal of this project is to extend the LightSABRE algorithm by integrating error-sensitive heuristics inspired by our previous work. The enhancements include:

- Incorporating error rates (CNOT gate errors, single-qubit gate errors, and readout errors) into routing decisions.
- Prioritizing critical paths and connectivity-aware edge selection.
- Adding dynamic calibration updates to account for real-time hardware noise.

## Roadmap

### 1. Modify `neighbor_table.rs`

#### What We Did:
- Replaced binary connectivity representation with error-aware weights using the formula: `-ln(1 - error_rate)`.
- Added functionality to fetch daily calibration data from backend properties and cache it locally.

#### Why We Did It:
- To penalize high-error connections during routing decisions, ensuring that paths with lower noise are prioritized.

### 2. Update `heuristic.rs`

#### What We Did:
- Implemented a new heuristic called **Connectivity-Aware Error Sensitive (CAES)**, which combines:
  - Connectivity information.
  - Error sensitivity.
  - Critical path prioritization (gates with high descendant counts are weighted higher).
- Updated scoring logic in `choose_best_swap` to compute a combined heuristic score.

#### Why We Did It:
- To balance connectivity constraints with fidelity optimization while ensuring that critical gates are routed through low-error paths.

### 3. Adjust `route.rs`

#### What We Did:
- Updated pathfinding algorithms (e.g., Dijkstra) to use error rates as edge weights instead of topological distances.
- Modified swap scoring to include penalties for frequently swapped qubits (using decay-based penalties).
- Incorporated readout errors into scoring mechanisms during routing trials.

#### Why We Did It:
- To minimize accumulated noise while avoiding overuse of specific physical qubits.

### 4. Enhance `layout.rs`

#### What We Did:
- Updated the initial layout generation process (`compute_dense_starting_layout`) to prioritize physical qubits with lower error rates during logical-to-physical mapping.

#### Why We Did It:
- Starting with a low-error initial layout reduces SWAP operations later and improves circuit fidelity from the beginning.

### 5. Refine `layer.rs`

#### What We Did:
- Penalized swaps involving high-readout-error qubits during scoring by adding readout error penalties in `FrontLayer::score()`.

#### Why We Did It:
- To improve measurement fidelity by avoiding high-readout-error qubits in critical operations.

### 6. Validation Tests

#### What We Did:
- Added unit tests to verify error-aware path selection, heuristic scoring, and readout error handling.
- Benchmarked the implementation using standard circuits like QFT and Grover's algorithm on real quantum hardware.

#### Why We Did It:
- To ensure that the modifications improve circuit fidelity without introducing runtime inefficiencies or errors in the mapping process.

## Key Features

### Dynamic Calibration Updates
- Error rates are fetched daily from backend properties, ensuring up-to-date noise profiles without runtime overhead.

### Critical Path Prioritization
- Gates along critical paths are prioritized in heuristic scoring, reducing overall circuit depth and accumulated errors.

### Hybrid Heuristics
- CAES heuristics are combined with decay-based penalties for frequently swapped qubits, balancing routing efficiency and fidelity.

## How It Works

1. The adjacency matrix in `neighbor_table.rs` is updated with error-aware weights based on real-time calibration data.
2. The CAES heuristic in `heuristic.rs` computes a combined score for each routing decision, balancing connectivity, error sensitivity, and critical path importance.
3. Pathfinding algorithms in `route.rs` use error-weighted distances to determine optimal SWAP operations.
4. Initial layouts in `layout.rs` prioritize low-error physical qubits for logical-to-physical mapping.
5. Scoring mechanisms in `layer.rs` penalize high-readout-error qubits during swap decisions.

## Validation Results

- Improved success probabilities compared to the original LightSABRE implementation.
- Reduced SWAP counts and circuit depths on standard benchmarks like QFT circuits.
- Demonstrated scalability and accuracy on real quantum hardware using IBM devices.

## How to Use

1. Clone this repository:
    ```
    git clone https://github.com/Legend-2727/qiskit.git
    cd qiskit
    ```

2. Build the Rust components:
    ```
    cargo build
    ```

3. Run unit tests:
    ```
    cargo test
    ```

4. Use the enhanced LightSABRE algorithm for transpiling quantum circuits with error-sensitive routing.
