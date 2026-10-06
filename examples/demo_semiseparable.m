%DEMO_SEMISEPARABLE  Test semiseparable block-encoding circuit
%
% This script generates random vectors u and v, builds the semiseparable
% matrix S = tril(u*v') + triu(v*u',1), constructs the corresponding
% quantum circuit, and compares the extracted block with the classical
% matrix.

clear;
clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

%% Problem size
N = 8;
rng(42);

%% Random normalized vectors
u = randn(N,1);
u = u / norm(u);

v = randn(N,1);
v = v / norm(v);

%% Classical semiseparable matrix
S = tril(u*v') + triu(v*u', 1);

%% Build quantum circuit
circ = build_semiseparable_circuit(u, v);

%% Visualize circuit
circ.draw;

%% Extract matrix from circuit
U = circ.matrix;

%% Extract encoded block
mat = U(1:N, 1:N) * 2 * sqrt(N);

%% Error analysis
err = norm(mat - S);

%% Display results
disp('Classical matrix S:');
disp(S);

disp('Extracted block from circuit:');
disp(mat);

disp('Error ||mat - S||:');
disp(err);