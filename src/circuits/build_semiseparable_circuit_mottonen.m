function circ = build_semiseparable_circuit_mottonen(v,u)
%BUILD_SEMISEPARABLE_CIRCUIT_MOTTONEN Same as BUILD_SEMISEPARABLE_CIRCUIT, but
% prepares u and v with QCLAB's InitializeStateVector (Mottonen state
% preparation) instead of a dense unitary. Constructs the block-encoding
% circuit for a
% semiseparable matrix.
%
% Given normalized vectors u and v, this routine constructs the circuit
% encoding the semiseparable matrix
%
%   S = tril(u*v') + triu(v*u', 1).
%
% Inputs:
%   u - Column vector of length N = 2^n.
%   v - Column vector of length N = 2^n.
%
% Output:
%   circ - qclab.QCircuit implementing the block encoding.

    %% Validate inputs

    N = length(u);
    n = log2(N);

    assert(mod(n,1) == 0, ...
        'The length of u and v must be a power of two.');
    assert(length(v) == N, ...
        'Vectors u and v must have the same length.');

    % Normalize 
    u = u / norm(u);
    v = v / norm(v);

    %% Allocate circuit

    numQubits = 3*n + 2;
    circ = qclab.QCircuit(numQubits);

    %% State preparation


circ.InitializeStateVector(v, 0);
circ.InitializeStateVector(u, n+1);

    %% State-to-operator transformation

    initCirc = init_circ(n);
    circ.push_back(initCirc);

    %% Z_g operator

    Z_sel = Sel_gate(n);
    
    circ.push_back(Z_sel);
    %circ.barrier(true);
    %% Z_s operator

    circ.push_back(Fix_gate(n));

    %% Final Hadamard

    circ.push_back(qclab.qgates.Hadamard(n));

    %% Pairwise CNOT layer

    circ.push_back(Flip_gate(n));
    %% Hadamard layer

    for i = 0:n-1
        circ.push_back(qclab.qgates.Hadamard(i));
    end

end