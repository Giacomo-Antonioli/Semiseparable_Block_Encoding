N=4;
n=log2(N);
%%%
x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);
S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix


X = sym('x', [N 1]);
Y = sym('y', [N 1]);

UX=repmat(sym('g'), N,N);
z=zeros(N,N);
UX=sym(z);
UX(:,1)=X;
UY=repmat(sym('g'), N,N);
UY=sym(z);
UY(:,1)=Y;
AB=kron(eye(2^(n+1)),kron(UY,kron(eye(2),UX)));
tot_wires = 3*n + 2;
circ = qclab.QCircuit(tot_wires);

for i=0:n-1
    circ.push_back(qclab.qgates.SWAP(i,2*n+2+i));
end
mat=circ.matrix;
res=AB*sym(mat);
circ.draw
tot_wires = 3*n + 2;
circ = qclab.QCircuit(tot_wires);
initCirc=init_circ(n);
%initCirc.asBlock("P");

circ.push_back(initCirc);
tmpcirc = qclab.QCircuit(tot_wires);
tmpcirc.push_back(initCirc);
init_circ_mat=tmpcirc.matrix;
tmpcirc.draw
tmp1=init_circ_mat*res;
%spy(res);
initmat=circ.matrix;
circBB=BB_gate(n);
circ.push_back(circBB);
tmpcirc = qclab.QCircuit(tot_wires);
tmpcirc.push_back(circBB);
bbmat=tmpcirc.matrix;
tmp2=bbmat*tmp1;
tmpcirc.draw
tmpcirc.toTex
% spy(tmp2)
tmpcirc = qclab.QCircuit(tot_wires);
for i=0:n-1
    circ.push_back(qclab.qgates.Hadamard(i));
    tmpcirc.push_back(qclab.qgates.Hadamard(i));
end
tmpcirc.draw
% 
% circ.draw
% 
% 
finalmat=circ.matrix;
ffmaat=finalmat*res;
spy(ffmaat(1:N,1:N))
