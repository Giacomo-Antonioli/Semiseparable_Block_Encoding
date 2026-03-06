

N=8;
n=log2(N);

x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);
% x=pi/4*(1:8)/norm(pi/4*(1:8));
% y=pi/4*(1:8)/norm(pi/4*(1:8));
% x=x';
% y=y';
S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tot_wires = 3*n + 2;
circ = qclab.QCircuit(tot_wires);
circBB=BB_gate(n);
circ.push_back(circBB);
for i=0:n-1
    circ.push_back(qclab.qgates.Hadamard(i));
end

mat=circ.matrix;
initcircuit=qclab.QCircuit(3*n+2);
Ux=make_circuit(x,n+1);
Ux.asBlock("Uy");
Uy=make_circuit(y,2*n+2);
Uy.asBlock("Uy");

initcircuit.push_back(Ux);

initcircuit.push_back(Uy);

circ=init_circ(n);
circ.draw
initcircuit.push_back(init_circ(n))
tmp=initcircuit.matrix;
% finalmat=mat*tmp;
finalmat=mat*tmp;
myS=zeros(N,N);
myS=full(finalmat(1:N,1:N)*2^(n/2+1));
norm(abs(myS)-abs(S))

S
myS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spy(tmp)
% 
% tot_wires = 3*n + 2;
% 
% initcircuit=qclab.QCircuit(tot_wires);
% Ux=make_circuit(x,n+1);
% Ux.asBlock("Uy");
% Uy=make_circuit(y,2*n+2);
% Uy.asBlock("Uy");
% 
% initcircuit.push_back(Ux);
% 
% initcircuit.push_back(Uy);
% 
% hswap=init_circ(n);
% hswap.asBlock("hswap")
% 
% initcircuit.push_back(hswap)
% 
% circ = qclab.QCircuit(tot_wires);
% circ.push_back(initcircuit)
% 
% circBB=BB_gate(n);
% circ.push_back(circBB);
% for i=0:n-1
%     circ.push_back(qclab.qgates.Hadamard(i));
% end
% 
% 
% circ.draw
% mat=circ.matrix;
% 
% 
% 
% finalmat=mat*tmp;
% myS=zeros(N,N);
% 
% 
