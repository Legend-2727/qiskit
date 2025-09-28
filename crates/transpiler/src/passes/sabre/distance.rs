// This code is part of Qiskit.
//
// (C) Copyright IBM 2025
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use fixedbitset::FixedBitSet;
use ndarray::{Array2, ArrayViewMut1, Axis};
use rayon_cond::CondIterator;
use rustworkx_core::petgraph::visit::{IntoNeighbors, NodeCompactIndexable};
use std::collections::HashMap;

// The implementation of `distance_matrix` was forked from Rustworkx at its commit 30f29079eeae,
// from the file `src/shortest_path/distance_matrix.rs` (as `compute_distance_matrix`). Its licence
// terms are:
//
//      Licensed under the Apache License, Version 2.0 (the "License"); you may
//      not use this file except in compliance with the License. You may obtain
//      a copy of the License at
//
//          http://www.apache.org/licenses/LICENSE-2.0
//
//      Unless required by applicable law or agreed to in writing, software
//      distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
//      WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
//      License for the specific language governing permissions and limitations
//      under the License.
//
// The implementation was modified to be generic over `petgraph::visit` traits, throw away the
// handling of `StableGraph`, minimise allocations in the `bfs_traversal` call, and use a bitmap
// data structure to track "seenness" rather than `HashSet`.

// This file may be obsoleted and/or upstreamed into Rustworkx in the future.

pub fn distance_matrix<G>(graph: G, parallel_threshold: usize, null_value: f64) -> Array2<f64>
where
    G: NodeCompactIndexable + IntoNeighbors,
{
    let n = graph.node_count();
    let neighbors = (0..n)
        .map(|index| {
            graph
                .neighbors(graph.from_index(index))
                .map(|neighbor| graph.to_index(neighbor))
                .collect::<FixedBitSet>()
        })
        .collect::<Vec<_>>();
    let bfs_traversal = |start: usize, mut row: ArrayViewMut1<f64>| {
        let mut distance = 0.0;
        let mut seen = FixedBitSet::with_capacity(n);
        let mut next = FixedBitSet::with_capacity(n);
        let mut cur = FixedBitSet::with_capacity(n);
        cur.put(start);
        while !cur.is_clear() {
            next.clear();
            for found in cur.ones() {
                row[[found]] = distance;
                next |= &neighbors[found];
            }
            seen.union_with(&cur);
            next.difference_with(&seen);
            distance += 1.0;
            ::std::mem::swap(&mut cur, &mut next);
        }
    };
    let mut out = Array2::from_elem((n, n), null_value);
    CondIterator::new(out.axis_iter_mut(Axis(0)), n >= parallel_threshold)
        .enumerate()
        .for_each(|(index, row)| bfs_traversal(index, row));
    out
}

/// Convert error probability to edge weight using CAES formula.
/// 
/// Args:
///     error_prob: Error probability (0.0 to 1.0, exclusive of 1.0)
/// 
/// Returns:
///     Edge weight calculated as -3.0 * ln(1.0 - clamp(p, 0 ≤ p < 1))
///     Represents cost of ~3 CNOTs for a SWAP operation
fn error_prob_to_weight(error_prob: f64) -> f64 {
    let clamped = error_prob.max(0.0).min(0.999999); // Avoid ln(0)
    -3.0 * (1.0 - clamped).ln()
}

/// Compute error-weighted all-pairs shortest paths using Floyd-Warshall algorithm.
/// This implements the CAES (Calibration-Aware Error-weighted SABRE) distance calculation.
/// 
/// Args:
///     graph: The coupling graph structure  
///     edge_errors: Map from (node1, node2) tuples to error rates (0.0 to 1.0)
/// 
/// Returns:
///     Distance matrix where distances are weighted by edge error rates.
///     Uses a large sentinel value (1e9) instead of infinity to avoid NaN issues.
pub fn compute_error_weighted_apsp<G>(
    graph: G,
    edge_errors: &HashMap<(usize, usize), f64>,
) -> Array2<f64>
where
    G: NodeCompactIndexable + IntoNeighbors,
{
    let n = graph.node_count();
    let large_value = 1e9; // Sentinel for infinity to avoid NaN/Inf issues
    let mut dist = Array2::from_elem((n, n), large_value);
    
    // Initialize diagonal to 0
    for i in 0..n {
        dist[[i, i]] = 0.0;
    }
    
    // Initialize direct edges with error weights
    for i in 0..n {
        let node_i = graph.from_index(i);
        for neighbor in graph.neighbors(node_i) {
            let j = graph.to_index(neighbor);
            let edge_key = if i < j { (i, j) } else { (j, i) };
            
            let weight = if let Some(&error_prob) = edge_errors.get(&edge_key) {
                error_prob_to_weight(error_prob)
            } else {
                // No calibration data - treat as unusable (high cost) 
                large_value
            };
            
            // Only update if we have a valid (non-infinite) weight
            if weight < large_value {
                dist[[i, j]] = weight;
            }
        }
    }
    
    // Floyd-Warshall all-pairs shortest paths
    for k in 0..n {
        for i in 0..n {
            for j in 0..n {
                if dist[[i, k]] < large_value && dist[[k, j]] < large_value {
                    let new_dist = dist[[i, k]] + dist[[k, j]];
                    if new_dist < dist[[i, j]] {
                        dist[[i, j]] = new_dist;
                    }
                }
            }
        }
    }
    
    dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use approx::assert_abs_diff_eq;
    use rustworkx_core::petgraph::{Graph, Undirected};

    /// Test the error probability to weight conversion formula
    #[test]
    fn test_error_prob_to_weight() {
        // Test p=0.01 → ~0.03015 (approximately 3 * 0.01005)
        let w = error_prob_to_weight(0.01);
        assert_abs_diff_eq!(w, 0.030151, epsilon = 1e-4);
        
        // Test p=0.0 → 0.0
        let w = error_prob_to_weight(0.0);
        assert_abs_diff_eq!(w, 0.0, epsilon = 1e-6);
        
        // Test p=0.5 → ~2.079
        let w = error_prob_to_weight(0.5);
        assert_abs_diff_eq!(w, 2.079442, epsilon = 1e-5);
    }

    /// Test APSP (All-Pairs Shortest Path) with error-weighted edges on a simple line graph
    #[test]
    fn test_apsp_error_weighted_line() {
        // Create a simple 3-node line graph: 0 -- 1 -- 2
        let edges = vec![(0u32, 1u32), (1u32, 2u32)];
        let graph = Graph::<(), (), Undirected, u32>::from_edges(&edges);
        
        // Define error rates: edge (0,1) has 1% error, edge (1,2) has 5% error
        let mut edge_errors = HashMap::new();
        edge_errors.insert((0, 1), 0.01);  // 1% error → weight ~0.0301
        edge_errors.insert((1, 2), 0.05);  // 5% error → weight ~0.1541
        
        let distance_matrix = compute_error_weighted_apsp(&graph, &edge_errors);
        
        // Verify distances
        // (0,0), (1,1), (2,2) should be 0
        assert_abs_diff_eq!(distance_matrix[[0, 0]], 0.0, epsilon = 1e-6);
        assert_abs_diff_eq!(distance_matrix[[1, 1]], 0.0, epsilon = 1e-6);
        assert_abs_diff_eq!(distance_matrix[[2, 2]], 0.0, epsilon = 1e-6);
        
        // (0,1) and (1,0) should be weight of error 0.01 → ~0.0301
        let expected_01 = error_prob_to_weight(0.01);
        assert_abs_diff_eq!(distance_matrix[[0, 1]], expected_01, epsilon = 1e-6);
        assert_abs_diff_eq!(distance_matrix[[1, 0]], expected_01, epsilon = 1e-6);
        
        // (1,2) and (2,1) should be weight of error 0.05 → ~0.1541
        let expected_12 = error_prob_to_weight(0.05);
        assert_abs_diff_eq!(distance_matrix[[1, 2]], expected_12, epsilon = 1e-6);
        assert_abs_diff_eq!(distance_matrix[[2, 1]], expected_12, epsilon = 1e-6);
        
        // (0,2) and (2,0) should be sum: ~0.0301 + ~0.1541 = ~0.1842
        let expected_02 = expected_01 + expected_12;
        assert_abs_diff_eq!(distance_matrix[[0, 2]], expected_02, epsilon = 1e-6);
        assert_abs_diff_eq!(distance_matrix[[2, 0]], expected_02, epsilon = 1e-6);
    }

    /// Test readout penalty calculation
    #[test]
    fn test_readout_penalty() {
        // Test with omega=1.0 and two qubits with readout errors 0.02 and 0.03
        let readout_a = 0.02;
        let readout_b = 0.03;  
        let omega = 1.0;
        
        let penalty = omega * (readout_a + readout_b);
        assert_abs_diff_eq!(penalty, 0.05, epsilon = 1e-6);
        
        // Test with different omega
        let omega = 2.5;
        let expected = 2.5 * (0.02 + 0.03);
        let actual = omega * (readout_a + readout_b);
        assert_abs_diff_eq!(actual, expected, epsilon = 1e-6);
    }

    /// Test fallback behavior when no error data is provided
    #[test] 
    fn test_fallback_no_errors() {
        // Create a simple 3-node line graph: 0 -- 1 -- 2  
        let edges = vec![(0u32, 1u32), (1u32, 2u32)];
        let graph = Graph::<(), (), Undirected, u32>::from_edges(&edges);
        
        // Empty error map should result in no valid paths (large values)
        let edge_errors = HashMap::new();
        let error_matrix = compute_error_weighted_apsp(&graph, &edge_errors);
        
        // All non-diagonal entries should be large_value since no error data means unusable edges
        let large_value = 1e9;
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    assert_abs_diff_eq!(error_matrix[[i, j]], 0.0, epsilon = 1e-6);
                } else {
                    assert_abs_diff_eq!(error_matrix[[i, j]], large_value, epsilon = 1e6);
                }
            }
        }
    }
}
