%TEST_CORRECTNESS_SEMISEPARABLE
% Verifies correctness of semiseparable block-encoding construction

clear;
clc;

%% Problem size
N = 4;
rng(42);

%% Generate normalized vectors (baseline test)
u = randn(N,1);
u = u / norm(u);

v = randn(N,1);
v = v / norm(v);

%% Classical target matrix
S = tril(u*v') + triu(v*u', 1);

%% Build block-encoding circuit
circ = build_semiseparable_circuit(u, v);

%% Extract full unitary
U = circ.matrix;

%% ============================================================
% (A) UNITARITY CHECK
%% ============================================================

I = eye(size(U));
unitarity_error = norm(U' * U - I, 2);

fprintf('\n=== Unitarity check ===\n');
fprintf('||U^†U - I||_2 = %.3e\n', unitarity_error);

%% ============================================================
% (B) BLOCK EXTRACTION CHECK
%% ============================================================

block_dim = N;
scale = 2 * sqrt(N);

U_block = U(1:block_dim, 1:block_dim);
S_tilde = scale * U_block;

block_error = norm(S_tilde - S, 2);

fprintf('\n=== Block encoding check ===\n');
fprintf('||S - S_tilde||_2 = %.3e\n', block_error);

%% ============================================================
% Summary
%% ============================================================

tol = 1e-10;

if unitarity_error < tol && block_error < tol
    fprintf('\nRESULT: PASS ✅\n');
else
    fprintf('\nRESULT: FAIL ❌\n');
end