import numpy as np
import pytest
from qiskit.quantum_info import Operator

import semiseparable_qiskit as s


@pytest.mark.parametrize("n", [1, 2, 3])
def test_mottonen_prepares_state(n):
    rng = np.random.default_rng(n)
    psi = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
    psi /= np.linalg.norm(psi)
    out = Operator(s.mottonen_state_prep_qiskit(psi)).data[:, 0]
    assert abs(np.vdot(psi, out)) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("n", [1, 2])
def test_full_circuit_is_unitary(n):
    circ = s.full_circuit(n)
    assert circ.num_qubits == 3 * n + 2
    U = Operator(circ).data
    assert np.linalg.norm(U.conj().T @ U - np.eye(U.shape[0])) < 1e-10


def test_demo_runs(capsys):
    s.demo(n=2, seed=0)
    last = capsys.readouterr().out.strip().splitlines()[-1]
    assert float(last) < 1e-10
