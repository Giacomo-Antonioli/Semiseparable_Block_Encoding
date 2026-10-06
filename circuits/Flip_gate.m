function circ = Flip_gate(n,asBlockFlag)
%FLIP_GATE  Apply the fixed CNOT ladder entangling layer.
%
% This gate implements the following pattern:
%
%   for i = 0:n-1
%       CNOT(i, n+i+1)
%   end
%
% It creates a structured entangling "flip" between two register blocks
% in a 3n+2 wire layout.
%
% INPUT
%   n
%       Number of qubits in the logical block.
%       Must be a positive integer (n > 0).
%
% OUTPUT
%   circ
%       qclab.QCircuit implementing the fixed entangling layer.
%
    if nargin < 2
        asBlockFlag = false;
    end
    % --- Assertions ---
    assert(n > 0, 'n must be greater than zero.');

    % total wires must match global encoding layout
    tot_wires = 3*n + 2;

    circ = qclab.QCircuit(tot_wires);

    for i = 0:n-1
        circ.push_back(qclab.qgates.CNOT(i, n + i + 1));
    end
  if asBlockFlag
    circ.asBlock("Z_{flip}");
  end

end