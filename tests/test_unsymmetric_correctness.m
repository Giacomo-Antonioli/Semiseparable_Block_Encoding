%TEST_UNSYMMETRIC_CORRECTNESS
% Verifies correctness of the unsymmetric semiseparable block-encoding
% construction S = tril(u*v') + triu(x*y',1).

clear;
clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();

%% Problem size
N = 8;
rng(42);

%% Random generators and classical target matrix
[S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);

%% Build block-encoding circuit with the normalized generators
circ = build_unsymmetric_semiseparable_circuit(u/nu, v/nv, x/nx, y/ny, ...
                                               nu, nv, nx, ny);

%% Extract the encoded block by simulation
dim = 2^circ.nbQubits;
psi = zeros(dim, 1);
S_tilde = zeros(N);
for j = 1:N
    psi(:) = 0;
    psi(j) = 1;
    sim = circ.simulate(psi);
    out = sim.states;
    S_tilde(:, j) = out(1:N);
end

%% Compare against S / alpha, alpha = sqrt(N)*(||u|| ||v|| + ||x|| ||y||)
alpha = sqrt(N)*(nu*nv + nx*ny);
block_error = norm(S/alpha - S_tilde, 2);

fprintf('\n=== Block encoding check ===\n');
fprintf('||S/alpha - S_tilde||_2 = %.3e\n', block_error);

tol = 1e-10;

if block_error < tol
    fprintf('\nRESULT: PASS\n');
else
    fprintf('\nRESULT: FAIL\n');
end
