function qc = state_prep_from_vector(v)
% STATE_PREP_FROM_VECTOR
% Prepare an n-qubit state from vector v (complex)
%
% Input:
%   v - column or row vector of length 2^n
% Output:
%   qc - qclab.QCircuit implementing the state preparation
    v = v(:);  % ensure column
    dim = length(v);
    n = log2(dim);
    assert(abs(n-round(n)) < 1e-12, 'Vector length must be 2^n');
    n = round(n);
    
    % Normalize input
    v = v / norm(v);
    
    % Initialize circuit
    qc = qclab.QCircuit(n);
    
    % Start recursive state preparation
    qc = recursive_prepare(qc, v, 0:(n-1));
end

%% Recursive helper function
function qc = recursive_prepare(qc, v, qubits)
    N = length(v);
    
    if N == 1
        return;
    end
    
    top = qubits(1);
    rest = qubits(2:end);
    n = N / 2;
    v0 = v(1:n);
    v1 = v(n+1:end);
    r0 = norm(v0);
    r1 = norm(v1);
    
    % RY rotation for top qubit
    if r0 + r1 > 1e-12
        theta = 2 * atan2(r1, r0);
        theta
        top
        qc.push_back(qclab.qgates.RotationY(int32(top), theta));
    end
    
    % Recurse on |0> branch
    if r0 > 1e-12 && ~isempty(rest)
        qc = recursive_prepare(qc, v0 / r0, rest);
    end
    
    % Recurse on |1> branch (controlled)
    if r1 > 1e-12 && ~isempty(rest)
        subqc = qclab.QCircuit(length(rest));
        subqc = recursive_prepare(subqc, v1 / r1, rest);
        
        % Apply subcircuit controlled on top qubit being |1>
        for g = 1:subqc.nGates
            gate = subqc.gates{g};
            qc.push_back(qclab.qgates.QControlledGate(gate, int32(top), 1));
        end
    end
    
    % Handle last qubit (rest empty) - apply phase correction
    if isempty(rest)
        for k = 1:N
            amp = v(k);
            phi = angle(amp);
            if abs(phi) > 1e-12
                qc.push_back(qclab.qgates.RotationZ(int32(top), phi));
            end
        end
    end
end

%% Test code
% Number of qubits
n = 3;

% Random complex vector
v = randn(2^n, 1) + 1i * randn(2^n, 1);
v = v / norm(v);

% Build circuit
qc = state_prep_from_vector(v);

% Get unitary matrix
U = qc.matrix();

% Apply to |0...0>
init = zeros(2^n, 1);
init(1) = 1;
statevec = U * init;

% Fidelity check
fidelity = abs(v' * statevec)^2;
disp(['Fidelity = ' num2str(fidelity)]);  % should be ~1

% Compare amplitudes
disp('Target vs output:')
disp([v statevec])