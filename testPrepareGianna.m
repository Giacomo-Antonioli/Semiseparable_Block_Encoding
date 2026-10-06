% --- verifica ---
n = 6; N = 2^n; theta = 3;
[circU, circV] = prepareUV(n, theta);

simU = circU.simulate( repmat('0',1,n) );  psiU = simU.states;
simV = circV.simulate( repmat('0',1,n) );  psiV = simV.states;

t = (0:N-1)'/N;
u_target = exp(theta*t);   u_target = u_target/norm(u_target);
v_target = exp(-theta*t);  v_target = v_target/norm(v_target);

fprintf('psiU - u_target = %.3e\n', norm(psiU(:) - u_target(:)));
fprintf('psiV - v_target = %.3e\n', norm(psiV(:) - v_target(:)));