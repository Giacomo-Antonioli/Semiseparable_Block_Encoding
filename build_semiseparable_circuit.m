function circ = build_semiseparable_circuit(u, v)
% BUILD_SEMISEPARABLE_CIRCUIT  Constructs the block-encoding circuit for a
% semiseparable matrix defined by vectors u and v.
%
% Inputs:
%   u - normalized column vector (length N, power of 2)
%   v - normalized column vector (length N, power of 2)
%
% Output:
%   circ - qclab.QCircuit encoding the semiseparable matrix S = tril(u*v') + triu(v*u', 1)

    % Infer size and number of qubits
    N = length(u);
    n = log2(N);

    assert(mod(n, 1) == 0, 'Length of u and v must be a power of 2.');
    assert(length(v) == N, 'u and v must have the same length.');

    % Normalize inputs (defensive copy)
    u = u / norm(u);
    v = v / norm(v);

    % Build circuit
    tot_wires = 3*n + 2;
    circ = qclab.QCircuit(tot_wires);

    % Initialize state vectors (note: u/v roles match original x/y mapping)
    Uv = qclab.InitializeStatevector(v, n+1);
    Uv.asBlock("SP(v)");

    Uu = qclab.InitializeStatevector(u);
    Uu.asBlock("SP(u)");

    circ.push_back(Uv);
    circ.push_back(Uu);

    % Append initialization circuit
    initCirc = init_circ(n);
    initCirc.asBlock("P")
    circ.push_back(initCirc);

    % Append BB gate
    circBB = BB_gate(n);
    circ.push_back(circBB);

    % Append Hadamard layer
    for i = 0:n-1
        circ.push_back(qclab.qgates.Hadamard(i));
    end

end