%SIDE_BY_SIDE_MOTTONEN_QCLAB
% Small, human-readable check: build the circuit for a chosen psi, print
% psi next to the first column of circ.matrix, side by side. The global
% phase is now baked into the circuit itself, so the two columns should
% match directly (no phase alignment needed).

clear;
clc;

%% ============================================================
% CHOOSE psi HERE
%% ============================================================

n = 3;               % number of qubits
N = 2^n;

rng(0);
psi = randn(N, 1) + 1i * randn(N, 1);
psi = psi / norm(psi);

% To test your own state instead, just set psi directly, e.g.:
% psi = [1; 1; 1; 1; 1; 1; 1; 1] / sqrt(8);

%% ============================================================
% BUILD CIRCUIT AND EXTRACT FIRST COLUMN
%% ============================================================

circuit = mottonen_state_prep_qclab(psi);
U = circuit.matrix;

prepared = U(:, 1);   % = circuit.matrix * e1

%% ============================================================
% SIDE-BY-SIDE PRINTOUT
%% ============================================================

fprintf('n = %d, N = %d\n\n', n, N);

fprintf('%4s | %22s | %22s | %10s\n', ...
    'idx', 'psi (input)', 'circ.matrix(:,1)', '|diff|');
fprintf('%s\n', repmat('-', 1, 68));

for i = 1:N
    fprintf('%4d | %9.4f %+8.4fi | %9.4f %+8.4fi | %10.2e\n', ...
        i-1, ...
        real(psi(i)), imag(psi(i)), ...
        real(prepared(i)), imag(prepared(i)), ...
        abs(prepared(i) - psi(i)));
end

fprintf('\nmax |circ.matrix(:,1) - psi| = %.3e\n', max(abs(prepared - psi)));