N=8;
n=log2(N);
rng(42)
%%%
%
% cercare di mettere il filo dei data in fondo, cosi' e' in cima
%%%
x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);

S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix
e0=[1;0]; e1=[0;1];
j=zeros(N, 1);
j(2)=1;
phi_init=  1/sqrt(2)*kron(e0, kron(kron(y, e0), x))+1/sqrt(2)*kron(e1, kron(kron(x, e0), y))



t=qclab.QCircuit(4)
t.push_back(qclab.qgates.SWAP(0,3))
c=qclab.QCircuit(2*n+2)
c.push_back(qclab.qgates.Hadamard(0))
c.push_back(qclab.qgates.MCMatrixGate(0,[1:4],t.matrix))
c.push_back(qclab.qgates.MCMatrixGate(0,[2:5],t.matrix))