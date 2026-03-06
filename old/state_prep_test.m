n = 3;
N=2^n;

% Build psi = (1 + [0..2^n-1]) + i(1 + [0..2^n-1])
%psi = (1 + (0:(2^n - 1))).' + 1i * (1 + (0:(2^n - 1))).';
psi = randn(N, 1); psi= psi/ norm(psi);
psi=psi';
disp('|ψ>');
disp(psi);

% Build the circuit (your function)
qc = make_circuit(psi);

mat=qc.matrix;
%%%%%%%%%%%%%TO MODIFY
% You must implement how to simulate the circuit in MATLAB.
% For now, st = psi (placeholder).
% st = psi;   % REMOVE THIS once a simulator is implemented
st=mat*e(0,n);
st=st'
% Compute error
err = norm(st - psi);
fprintf('Error = %.12f\n', err);
% qc.draw
% qc.toTex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
