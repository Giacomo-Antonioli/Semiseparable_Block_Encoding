%% test_state_prep.m
% Test script for state preparation with custom MCPhase using utils.m

clear; clc;

% ------------------ Parameters ------------------
n = 3; % number of qubits
dim = 2^n;

% Generate a random normalized complex vector
v = randn(dim,1) + 1i*randn(dim,1);
v = v / norm(v);

disp('Target statevector:');
disp(v);

% ------------------ Build Circuit ------------------
% utils.m should include:
%   state_prep_from_vector_MCPhase
qc = state_prep_from_vector_MCPhase(v);

% ------------------ Simulate Circuit ------------------
statevec = qc.statevector();

disp('Circuit output statevector:');
disp(statevec);

% ------------------ Fidelity / Error ------------------
% Fidelity = |<target|output>|^2
fidelity = abs(v' * statevec)^2;
fprintf('Fidelity = %.16f\n', fidelity);

% Norm of difference as additional check
error = norm(statevec - v);
fprintf('Norm of difference = %.3e\n', error);

if error < 1e-12
    disp('SUCCESS: Circuit reproduces the target state.');
else
    disp('WARNING: Circuit does not match the target state.');
end

% ------------------ Optional: Display probabilities ------------------
probs = abs(statevec).^2;
disp('Output probabilities:');
disp(probs);

