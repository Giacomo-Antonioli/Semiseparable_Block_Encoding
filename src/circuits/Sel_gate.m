function circ = Sel_gate(n,asBlockFlag)
%SEL_GATE Construct the Z_sel operator.
%
% This gate implements a structured rightshift pattern used in the
% block-encoding construction.
%
% INPUT
%   n
%       Number of logical qubits (must be > 0)
%
% OUTPUT
%   circ
%       qclab.QCircuit implementing the Z_sel operator
%
% CIRCUIT STRUCTURE
%   The operator acts on a register of size:
%       2n + 2 wires
%
%   It applies a sequence of rightshift operations:
%       control wires: tot_wires-n-1 : tot_wires-n+i-1
%       shift index  : i
%       shift size   : 1
%
%   The pattern builds a forward propagation structure across the
%   upper register block.

    % --- Assertions ---
    assert(n > 0, 'n must be greater than zero.');
    if nargin < 2
        asBlockFlag = false;
    end
    % --- Build circuit ---
    tot_wires = 2*n + 2;
    circ = qclab.QCircuit(tot_wires, n);

    % --- Forward propagation structure ---


    for i = 0:n
        circ.push_back( ...
            rightshift( ...
                tot_wires, ...
                tot_wires-n-1 : tot_wires-n+i-1, ...
                i, ...
                1));
    end

    
    if asBlockFlag
        circ.asBlock("Z_{sel}");
    end
end