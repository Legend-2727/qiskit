// This code is part of Qiskit.
//
// (C) Copyright IBM 2022
//
// This code is licensed under the Apache License, Version 2.0. You may
// obtain a copy of this license in the LICENSE.txt file in the root directory
// of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
//
// Any modifications or derivative works of this code must retain this
// copyright notice, and modified files need to carry a notice indicating
// that they have been altered from the originals.

use crate::getenv_use_multiple_threads;
use ndarray::prelude::*;
use numpy::{PyArray2, PyArrayMethods, PyReadonlyArray2, ToPyArray};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use rayon::prelude::*;
use rustworkx_core::petgraph::prelude::*;
use smallvec::SmallVec;

use crate::nlayout::PhysicalQubit;

/// A simple container that contains a vector of vectors representing
/// neighbors of each node in the coupling map. Optionally stores weights,
/// distances, and readout errors if you want to incorporate hardware error rates.
#[pyclass(module = "qiskit._accelerate.sabre")]
#[derive(Clone, Debug)]
pub struct NeighborTable {
    // Each entry is a list of physical neighbors
    neighbors: Vec<SmallVec<[PhysicalQubit; 4]>>,
    // Optional example arrays for error-based logic:
    weights: Option<Array2<f64>>,
    distances: Option<Array2<f64>>,
    readout_errors: Option<Vec<f64>>,
}

impl NeighborTable {
    /// Rebuilds a Rust DiGraph from the table:
    pub fn coupling_graph(&self) -> DiGraph<(), ()> {
        let mut graph = DiGraph::new();
        // We assume contiguous node indices 0..neighbors.len()
        // so we can do a direct add of the correct number of nodes:
        let n = self.neighbors.len();
        for _ in 0..n {
            graph.add_node(());
        }
        for (u, nbrs) in self.neighbors.iter().enumerate() {
            for &phys in nbrs {
                let v = phys.index() as usize;
                graph.add_edge(NodeIndex::new(u), NodeIndex::new(v), ());
            }
        }
        graph
    }

    pub fn num_qubits(&self) -> usize {
        self.neighbors.len()
    }

    /// Example helper if you want to retrieve a distance:
    pub fn get_error_distance(&self, from: usize, to: usize) -> f64 {
        // Return 1.0 if there's no distance array or out-of-bounds
        if let Some(dist) = &self.distances {
            if from < dist.nrows() && to < dist.ncols() {
                return dist[[from, to]];
            }
        }
        1.0
    }

    pub fn get_readout_error(&self, qubit: usize) -> f64 {
        if let Some(errors) = &self.readout_errors {
            if qubit < errors.len() {
                return errors[qubit];
            }
        }
        0.0
    }
}

impl std::ops::Index<PhysicalQubit> for NeighborTable {
    type Output = [PhysicalQubit];

    fn index(&self, index: PhysicalQubit) -> &Self::Output {
        &self.neighbors[index.index()]
    }
}

#[pymethods]
impl NeighborTable {
    /// Construct a new neighbor table from adjacency_matrix, plus optional error map.
    #[new]
    #[pyo3(signature = (adjacency_matrix, error_map=None))]
    pub fn new(
        adjacency_matrix: Option<PyReadonlyArray2<f64>>,
        error_map: Option<Py<PyAny>>, // or Option<Bound<'_, PyAny>>
    ) -> PyResult<Self> {
        let run_in_parallel = getenv_use_multiple_threads();

        // Step 1: build neighbors
        let neighbors = if let Some(ref adj) = adjacency_matrix {
            let mat = adj.as_array();
            // Build the list of neighbors row by row
            if run_in_parallel {
                let rows: Vec<_> = mat.axis_iter(Axis(0)).collect();
                // mat.axis_iter(Axis(0))
                rows.into_par_iter()
                    .enumerate()
                    .map(|(_row_i, row)| {
                        let mut list = SmallVec::<[PhysicalQubit; 4]>::new();
                        for (col_i, &val) in row.iter().enumerate() {
                            if val != 0.0 {
                                list.push(PhysicalQubit::new(col_i as u32));
                            }
                        }
                        Ok::<_, PyErr>(list)
                    })
                    .collect::<Result<Vec<_>, _>>()?
            } else {
                let mut all_rows = Vec::with_capacity(mat.nrows());
                for row in mat.axis_iter(Axis(0)) {
                    let mut list = SmallVec::<[PhysicalQubit; 4]>::new();
                    for (col_i, &val) in row.iter().enumerate() {
                        if val != 0.0 {
                            list.push(PhysicalQubit::new(col_i as u32));
                        }
                    }
                    all_rows.push(list);
                }
                all_rows
            }
        } else {
            Vec::new()
        };

        // Step 2: optionally parse error_map, set up weights/distances/readout
        let mut weights: Option<Array2<f64>> = None;
        let mut distances: Option<Array2<f64>> = None;
        let mut readout_errors: Option<Vec<f64>> = None;

        if let Some(pyobj) = error_map {
            println!("[NeighborTable] Received error_map from Python!");
            // Prepare local variables
            let n = neighbors.len();
            let mut w = Array2::<f64>::from_elem((n, n), 1.0);
            let mut d = Array2::<f64>::from_elem((n, n), 1.0);
            let mut r = vec![0.0; n];

            // Access the Python dict inside the GIL closure
            Python::with_gil(|py| -> PyResult<()> {
                let bound_any = pyobj.into_bound(py);
                // Step 2: Downcast Bound<'py, PyAny> → Bound<'py, PyDict>
                let dict: Bound<'_, PyDict> = bound_any.downcast_into::<PyDict>()?;
                // Parse edges or readout errors
                for (k, v) in dict.iter() {
                    println!("[NeighborTable] Parsing key: {:?}, value: {:?}", k, v);
                    if let Ok((q1, q2)) = k.extract::<(usize, usize)>() {
                        // two-qubit error
                        let errval = v.extract::<f64>()?;
                        if q1 < n && q2 < n {
                            w[[q1, q2]] = -3.0 * (1.0 - errval).ln();
                            d[[q1, q2]] = w[[q1, q2]];
                        }
                    } else if let Ok(q) = k.extract::<usize>() {
                        // single-qubit readout error
                        let ro_err = v.extract::<f64>()?;
                        if q < n {
                            r[q] = ro_err;
                        }
                    }
                }
                println!("[NeighborTable] Done parsing error_map. w = {:?}, d = {:?}, r = {:?}", w, d, r);
                Ok(())
            })?;

            // Floyd–Warshall to fill out distances
            for k in 0..n {
                for i in 0..n {
                    for j in 0..n {
                        let alt = d[[i, k]] + d[[k, j]];
                        if alt < d[[i, j]] {
                            d[[i, j]] = alt;
                        }
                    }
                }
            }
            println!("[NeighborTable] Floyd–Warshall complete: distances = {:?}", d);
            println!("[NeighborTable] readout_errors = {:?}", r);
            weights = Some(w);
            distances = Some(d);
            readout_errors = Some(r);

        }

        Ok(Self {
            neighbors,
            weights,
            distances,
            readout_errors,
        })
    }

    /// Example: check if we have error data
    #[getter]
    pub fn has_error_data(&self) -> bool {
        self.weights.is_some() && self.distances.is_some() && self.readout_errors.is_some()
    }

    /// Example get distance
    #[pyo3(signature=(src, dst))]
    fn error_distance(&self, src: usize, dst: usize) -> f64 {
        self.get_error_distance(src, dst)
    }

    /// Example get readout
    #[pyo3(signature=(q))]
    fn qubit_readout_error(&self, q: usize) -> f64 {
        self.get_readout_error(q)
    }

    /// Return the picklable python state
    fn __getstate__(&self, py: Python<'_>) -> PyResult<Py<PyList>> {
        // We'll build up a PyList of length 4: [neighbors, weights, distances, readout]
        let state_list = PyList::empty(py);

        // 1) neighbors => convert each row to a python list of qubit indices
        let nlist = PyList::empty(py);
        for row in &self.neighbors {
            // convert smallvec to list of ints
            let row_ints: Vec<_> = row.iter().map(|q| q.index()).collect();
            let row_py = PyList::new(py, row_ints);
            let row_py = row_py?;
            nlist.append(row_py)?;
        }
        state_list.append(nlist)?;

        // 2) weights => if Some, store NxN as python array, else None
        if let Some(ref w) = self.weights {
            // convert to PyArray
            let pyarr = w.to_pyarray(py);
            state_list.append(pyarr)?;
        } else {
            state_list.append(py.None())?;
        }

        // 3) distances => if Some, store NxN, else None
        if let Some(ref d) = self.distances {
            let pyarr = d.to_pyarray(py);
            state_list.append(pyarr)?;
        } else {
            state_list.append(py.None())?;
        }

        // 4) readout => if Some, store python list
        if let Some(ref r) = self.readout_errors {
            let r_py = PyList::new(py, r)?;
            state_list.append(r_py)?;
        } else {
            state_list.append(py.None())?;
        }

        Ok(state_list.into())
    }

    /// Unpickle
    fn __setstate__(&mut self, state: Py<PyAny>) -> PyResult<()> {
        Python::with_gil(|py| {
            // First, store the Bound<'py, PyAny> in a local binding
            let bound_any = state.into_bound(py);
            // Then downcast that binding to a PyList
            let st_list = bound_any.downcast::<PyList>()?;

            if st_list.len() != 4 {
                return Err(PyValueError::new_err(
                    "Invalid state length in NeighborTable",
                ));
            }

            // neighbors
            let item_0 = st_list.get_item(0)?;
            let nlist = item_0.downcast::<PyList>()?;
            let mut new_neighbors = Vec::with_capacity(nlist.len());
            for sub in nlist.iter() {
                let sublist = sub.downcast::<PyList>()?;
                let mut rowvec = SmallVec::<[PhysicalQubit; 4]>::new();
                for item in sublist.iter() {
                    let i = item.extract::<u32>()?;
                    rowvec.push(PhysicalQubit::new(i));
                }
                new_neighbors.push(rowvec);
            }
            self.neighbors = new_neighbors;

            // weights
            let w_obj = st_list.get_item(1)?;
            if w_obj.is_none() {
                self.weights = None;
            } else {
                let arr = w_obj.downcast::<PyArray2<f64>>()?;
                self.weights = Some(arr.readonly().as_array().to_owned());
            }

            // distances
            let d_obj = st_list.get_item(2)?;
            if d_obj.is_none() {
                self.distances = None;
            } else {
                let arr = d_obj.downcast::<PyArray2<f64>>()?;
                self.distances = Some(arr.readonly().as_array().to_owned());
            }

            // readouts
            let r_obj = st_list.get_item(3)?;
            if r_obj.is_none() {
                self.readout_errors = None;
            } else {
                let r_py = r_obj.downcast::<PyList>()?;
                let mut vals = Vec::with_capacity(r_py.len());
                for item in r_py.iter() {
                    vals.push(item.extract::<f64>()?);
                }
                self.readout_errors = Some(vals);
            }

            Ok(())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use numpy::PyArray2;
    use pyo3::types::PyDict;
    // Test 1: Basic adjacency -> neighbor list
    #[test]
    fn test_basic_neighbors() {
        Python::with_gil(|py| {
            // Simple chain of 3 qubits: 0--1--2
            let mut adj_matrix = Array2::<f64>::zeros((3, 3));
            adj_matrix[[0, 1]] = 1.0;
            adj_matrix[[1, 0]] = 1.0;
            adj_matrix[[1, 2]] = 1.0;
            adj_matrix[[2, 1]] = 1.0;

            // Convert to PyReadonlyArray2
            let py_adj = PyArray2::from_array(py, &adj_matrix).readonly();
            // No error map
            let table = NeighborTable::new(Some(py_adj), None).unwrap();

            assert_eq!(table.num_qubits(), 3);
            // Qubit 0 should have neighbor 1
            assert_eq!(table.neighbors[0].len(), 1);
            assert_eq!(table.neighbors[0][0].index(), 1);
            // Qubit 1 should have neighbors 0, 2
            assert_eq!(table.neighbors[1].len(), 2);
        });
    }

    // Test 2: Check has_error_data == false and default fallback
    #[test]
    fn test_no_error_data() {
        Python::with_gil(|py| {
            let mut adj_matrix = Array2::<f64>::zeros((2, 2));
            adj_matrix[[0, 1]] = 1.0;
            adj_matrix[[1, 0]] = 1.0;

            let py_adj = PyArray2::from_array(py, &adj_matrix).readonly();
            // Construct with no error map
            let table = NeighborTable::new(Some(py_adj), None).unwrap();

            // Expect no error data
            assert!(!table.has_error_data());
            // Distances fallback to 1.0 for adjacent
            assert_eq!(table.get_error_distance(0, 1), 1.0);
            // Readout error fallback is 0.0
            assert_eq!(table.get_readout_error(0), 0.0);
        });
    }

    // Test 3: Partial error data => some edges/qubits are missing
    #[test]
    fn test_partial_error_data() {
        Python::with_gil(|py| {
            let mut adj_matrix = Array2::<f64>::zeros((3, 3));
            // 0--1--2 chain
            adj_matrix[[0, 1]] = 1.0;
            adj_matrix[[1, 0]] = 1.0;
            adj_matrix[[1, 2]] = 1.0;
            adj_matrix[[2, 1]] = 1.0;

            let error_map = PyDict::new(py);
            // Provide error only for edge (0,1)
            error_map.set_item((0, 1), 0.01).unwrap();
            // Omit (1,2) and readout errors => default to 0

            let py_adj = PyArray2::from_array(py, &adj_matrix).readonly();
            let table = NeighborTable::new(Some(py_adj), Some(error_map.into())).unwrap();

            // The path 0 -> 1 -> 2 includes:
            //   (0,1) => error = 0.01 => weight? = -3.0 * ln(1 - 0.01) ...
            //   (1,2) => error = 0.0 => weight = -3.0 * ln(1.0) => 0.0
            let dist_0_2 = table.get_error_distance(0, 2);
            // Make sure it's > 0 because of the (0,1) edge
            assert!(dist_0_2 > 0.0);
        });
    }

    // Test 4: Full error data => verifying readout errors and Floyd–Warshall
    #[test]
    fn test_error_weighted_distances() {
        Python::with_gil(|py| {
            // 3-qubit chain: 0--1--2
            let mut adj_matrix = Array2::<f64>::zeros((3, 3));
            adj_matrix[[0, 1]] = 1.0;
            adj_matrix[[1, 0]] = 1.0;
            adj_matrix[[1, 2]] = 1.0;
            adj_matrix[[2, 1]] = 1.0;

            let error_map = PyDict::new(py);
            // Edges
            error_map.set_item((0, 1), 0.01).unwrap(); // low error
            error_map.set_item((1, 2), 0.1).unwrap(); // higher error
                                                      // Readout
            error_map.set_item(0usize, 0.02).unwrap();
            error_map.set_item(1usize, 0.03).unwrap();
            error_map.set_item(2usize, 0.04).unwrap();

            let py_adj = PyArray2::from_array(py, &adj_matrix).readonly();
            let table = NeighborTable::new(Some(py_adj), Some(error_map.into())).unwrap();

            // Should have error data
            assert!(table.has_error_data());

            // Check readout
            assert!((table.get_readout_error(0) - 0.02).abs() < 1e-9);
            assert!((table.get_readout_error(1) - 0.03).abs() < 1e-9);
            assert!((table.get_readout_error(2) - 0.04).abs() < 1e-9);

            // Check distance 0->2 => should be sum of (0->1) + (1->2)
            let d_0_1 = table.get_error_distance(0, 1);
            let d_1_2 = table.get_error_distance(1, 2);
            let d_0_2 = table.get_error_distance(0, 2);

            // If you run Floyd–Warshall, it should find that path is 0->1->2
            // so distance(0,2) ~ d_0_1 + d_1_2
            let diff = (d_0_2 - (d_0_1 + d_1_2)).abs();
            assert!(
                diff < 1e-6,
                "0->2 should equal (0->1)+(1->2). Got {}",
                d_0_2
            );
        });
    }

    // Test 5: Check pickle/unpickle roundtrip (you already have a version, but leaving here)
    #[test]
    fn test_pickle_roundtrip() {
        Python::with_gil(|py| {
            // 1) Construct simple adjacency matrix
            let mut adj = Array2::<f64>::zeros((2, 2));
            adj[[0, 1]] = 1.0;
            adj[[1, 0]] = 1.0;
            let py_adj = PyArray2::from_array(py, &adj).readonly();

            // 2) Create original NeighborTable
            let table = NeighborTable::new(Some(py_adj), None).unwrap();

            // 3) Get the picklable state as Py<PyList> (automatically a Py<PyAny>)
            let state: Py<PyAny> = table.__getstate__(py).unwrap().into();

            // 4) Construct an empty NeighborTable
            let mut new_table = NeighborTable {
                neighbors: vec![],
                weights: None,
                distances: None,
                readout_errors: None,
            };

            // 5) Call __setstate__ with Py<PyAny>
            new_table.__setstate__(state).unwrap();

            // 6) Validate roundtrip worked
            assert_eq!(new_table.num_qubits(), 2);
            assert_eq!(new_table.neighbors[0].len(), 1);
        });
    }
}
