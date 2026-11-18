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
%   For i = 0:n-1 it applies a leftshift
%       target wires :  tot_wires-n : tot_wires-1-i
%       control wire :  n-1-i  (closed, state 1)
%
%   The placeholder wire n is not used as a control: both halves of the
%   shifted subvectors undergo the same shift and only the first half
%   enters the block encoding (see the paper, Section 4). An earlier
%   version also used wire n as an open control (state 0) in every shift;
%   the encoded top-left block is identical, only the part of the unitary
%   outside the block and the gate counts differ.
%
%   See also LEFTSHIFT.
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
            n-1-i, ...
            1);

        circ.push_back(tmp);
    end
    if asBlockFlag 
    
    circ.asBlock("Z_{fix}");
    end
end