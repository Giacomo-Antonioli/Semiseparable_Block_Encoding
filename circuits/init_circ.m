
function circ = init_circ(n)
  % init_circ: build a QCircuit with a SWAP-constructed matrix gate and Hadamard
  % Inputs:
  %   n - integer, number of control steps
  %
  % Output:
  %   circ - qclab.QCircuit of size (2*n + 2) with n logical qubits


  % Main circuit with (2*n + 2) physical qubits, logical dimension n
  circ = qclab.QCircuit(3 * n + 2);
% 
%   % Apply Hadamard on qubit 0
for i=0:n-1
    circ.push_back(qclab.qgates.SWAP(i,2*n+2+i));
end

 %circ.push_back(qclab.qgates.Hadamard(n));

  %Add MCMatrixGate for i = 0 .. n-1
 for i = 0 : n - 1
   circ.push_back(qclab.qgates.CSWAP(n,n+1+i,2*n+2+i));

 end


