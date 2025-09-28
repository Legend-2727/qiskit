from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag
from qiskit.transpiler import CouplingMap
from qiskit.transpiler.passes.routing import SabreSwap

def ops_list(dag):
    return [(n.op.name, tuple(dag.qubits.index(q) for q in n.qargs))
            for n in dag.topological_op_nodes()]

def test_parity_without_calibration():
    qc = QuantumCircuit(4)
    qc.cx(0, 3); qc.cx(1, 2)
    dag = circuit_to_dag(qc)
    cm = CouplingMap.from_line(4)

    base = SabreSwap(cm, heuristic='basic', seed=123).run(dag)
    err  = SabreSwap(cm, heuristic='error_aware', lambda_=0.6, omega=1.0, seed=123).run(dag)

    assert ops_list(base) == ops_list(err)
    assert base.depth() == err.depth()