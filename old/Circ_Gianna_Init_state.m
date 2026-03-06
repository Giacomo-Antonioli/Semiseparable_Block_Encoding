N=8;
n=log2(N);
x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);
S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix
e0=[1;0]; e1=[0;1];
j=zeros(N, 1);
j(2)=1;
phi_init=  1/sqrt(2)*kron(e0, kron(kron(y, e0), x))+1/sqrt(2)*kron(e1, kron(kron(x, e0), y));
t=qclab.QCircuit(n+2)
t.push_back(qclab.qgates.SWAP(0,n+1))
c=qclab.QCircuit(2*n+2,n)
c.push_back(qclab.qgates.Hadamard(0))
for i=0:n-1
c.push_back(qclab.qgates.MCMatrixGate(0,[1+i:n+2+i],t.matrix))
end

c.draw
init=kron([1;0],kron(y,kron([1;0],x)));
myres=c.matrix*init;
norm(abs(myres)-abs(phi_init))