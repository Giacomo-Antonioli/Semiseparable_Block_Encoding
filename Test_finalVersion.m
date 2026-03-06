N = 4; rng(42);
x = randn(N,1); x = x/norm(x);
y = randn(N,1); y = y/norm(y);
S = tril(x*y') + triu(y*x', 1);

circ = build_semiseparable_circuit(x, y);
n = log2(N);
finalmat = circ.matrix;
myS = full(finalmat(1:N, 1:N) * 2^(n/2+1));
circ.draw
fprintf("BE norm: %.6e\n", norm(abs(myS) - abs(S)));