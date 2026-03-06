N=8;
n=log2(N);
%%%
x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);
S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix

tot_wires = 3*n + 2;
circ = qclab.QCircuit(tot_wires);


for i=0:n-1
    circ.push_back(qclab.qgates.SWAP(i,2*n+2+i));
end
% Ux=make_circuit(y,n+1);
% Ux.asBlock("Uy");
% Uy=make_circuit(x,2*n+2);
% Uy.asBlock("Uy");
% 
% circ.push_back(Ux);
% 
% circ.push_back(Uy);
% 
% circ.push_back(init_circ(n));


initmat=circ.matrix;
circBB=BB_gate(n);
circ.push_back(circBB);



for i=0:n-1
    circ.push_back(qclab.qgates.Hadamard(i));
end



finalmat=circ.matrix;


myS=full(finalmat(1:N,1:N)*2^(n/2+1));
fprintf("BE norm")
norm(abs(myS)-abs(S))
