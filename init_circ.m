
function circ = init_circ(n)
  % init_circ: build a QCircuit with a SWAP-constructed matrix gate and Hadamard
  % Inputs:
  %   n - integer, number of control steps
  %
  % Output:
  %   circ - qclab.QCircuit of size (2*n + 2) with n logical qubits

  % Temporary circuit on (n + 2) qubits
  t = qclab.QCircuit(n + 2);

  % Build a SWAP between qubits 0 and (n + 1)
  t.push_back(qclab.qgates.SWAP(0, n + 1));

  % Main circuit with (2*n + 2) physical qubits, logical dimension n
  circ = qclab.QCircuit(2 * n + 2, n);

  % Apply Hadamard on qubit 0
  circ.push_back(qclab.qgates.Hadamard(0));

  % Add MCMatrixGate for i = 0 .. n-1
  for i = 0 : n - 1
    controls = 1 + i : n + 2 + i;  % vector of control indices
    circ.push_back(qclab.qgates.MCMatrixGate(0, controls, t.matrix,'CSWAP'));
  end
