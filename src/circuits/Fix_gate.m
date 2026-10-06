function circ = Fix_gate(n,asBlockFlag)
%FIX_GATE Construct the fixed shift-based entangling operator.
%
% This gate implements a structured leftshift pattern used in the
% block-encoding construction.
%
% INPUT
%   n
%       Number of logical qubits (must be > 0)
%
% OUTPUT
%   circ
%       qclab.QCircuit implementing the fixed shift operator
%
% CIRCUIT STRUCTURE
%   The operator acts on a register of size:
%       2n + 1 wires
%
%   It applies a sequence of leftshift operations:
%       control wires:  tot_wires-n : tot_wires-1-i
%       target wires :  [n-1-i, n]
%
    if nargin < 2
        asBlockFlag = false;
    end
    % --- Assertions ---
    assert(n > 0, 'n must be greater than zero.');

    % --- Build circuit ---
    tot_wires = 2*n + 1;

    circ = qclab.QCircuit(tot_wires, n+1);

    % --- Shift structure ---
    for i = 0:n-1

        tmp = leftshift( ...
            tot_wires, ...
            tot_wires - n : tot_wires - 1 - i, ...
            [n-1-i, n], ...
            [1, 0]);

        circ.push_back(tmp);
    end
    if asBlockFlag 
    
    circ.asBlock("Z_{fix}");
    end
end