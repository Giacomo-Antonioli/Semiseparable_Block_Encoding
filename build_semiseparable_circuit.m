function circ = build_semiseparable_circuit(v,u,asBlockFlag)
%BUILD_SEMISEPARABLE_CIRCUIT Construct the block-encoding circuit for a
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
    if nargin < 3
        asBlockFlag = false;
    end
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

    Uv = state_prep(v,0);
    if asBlockFlag
    Uv.asBlock("SP(v)");
    end
    Uu = state_prep(u,n+1);
    if asBlockFlag
    Uu.asBlock("SP(u)");
    end
    circ.push_back(qclab.qgates.Hadamard(n))
    circ.push_back(Uv);
    circ.push_back(Uu);

    %% State-to-operator transformation

    initCirc = init_circ(n);
    circ.push_back(initCirc);
    circ.barrier(true)
    %% Z_g operator

    Z_sel = Sel_gate(n,asBlockFlag);
    %Z_sel.barrier(true);

    circ.push_back(Z_sel);
    %circ.barrier(true);
    %% Z_s operator
circ.barrier(true)
    circ.push_back(Fix_gate(n,asBlockFlag));

    %% Final Hadamard

    circ.push_back(qclab.qgates.Hadamard(n));

    %% Pairwise CNOT layer

    circ.push_back(Flip_gate(n,asBlockFlag));
    %% Hadamard layer

    for i = 0:n-1
        circ.push_back(qclab.qgates.Hadamard(i));
    end

end