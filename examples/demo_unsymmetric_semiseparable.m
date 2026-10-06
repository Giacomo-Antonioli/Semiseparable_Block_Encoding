%DEMO_UNSYMMETRIC_SEMISEPARABLE  Test unsymmetric semiseparable block-encoding circuit
%


% clear;
% clc;
digits(64)
%% Problem size
N = 8;
addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
setup_paths();
rng(42);

[S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N);

% inv(S)
% Build quantum circuit
u=u/nu;
v=v/nv;
x=x/nx;
y=y/ny;
[circ,t1,t2] = build_unsymmetric_semiseparable_circuit(u, v, x, y,nu,nv,nx,ny);
dim = 2^circ.nbQubits;
psi = zeros(dim, 1);
S_tilde = zeros(N);

% 3. Extract matrix column by column via simulation
for j = 1:N
    j
    psi(:) = 0;
    psi(j) = 1;
    sim = circ.simulate(psi);
    out = sim.states;
    S_tilde(:, j) =   out(1:N);% *(nu*nv*nx*ny);/(t1*t2)sqrt(N) *
end


% 4. Compute absolute spectral error
err = norm(S/(sqrt(N)*(nu*nv+nx*ny))-S_tilde, 2);

%% Display results
% disp('Classical matrix S:');
disp(S/(sqrt(N)*(nu*nv+nx*ny)));

% disp('Extracted block from circuit:');
disp(S_tilde);

disp('Error ||S_tilde - S||:');
disp(err);
(S/(sqrt(N)*(nu*nv+nx*ny)))./S_tilde