N=8;
n=log2(N);
%%%
x = randn(N, 1); x = x / norm(x);
y = randn(N, 1); y = y / norm(y);
S = tril(x*y')+ triu(y*x', 1); %semiseparable matrix
e0=[1;0]; e1=[0;1];
j=zeros(N, 1);
j(2)=1;
phi_init=  1/sqrt(2)*kron(e0, kron(kron(y, e0), x))+1/sqrt(2)*kron(e1, kron(kron(x, e0), y));

N2=2*N;
Z2N=zeros(N2); Z2N=diag(ones(N2-1, 1), 1); Z2N(N2, 1)=1; %upshift
Zb=eye(N*N2);
for i=1:N-1
    Zb(i*N2+1:(i+1)*N2, i*N2+1:(i+1)*N2)=Z2N^i;
end
Zbi=eye(N*N2);
for i=0:N-1
    Zbi(i*N2+1:(i+1)*N2, i*N2+1:(i+1)*N2)=(Z2N')^(N-i);
end
spy(Zb)

Zg= blkdiag(Zb, Zbi);
ZN=zeros(N); ZN=diag(ones(N-1, 1), -1); ZN(1, N)=1; % small downshift
Zs=eye(N*N2);
for i=1:N-1
    Zs(i*N2+1:(i)*N2+N, i*N2+1:(i)*N2+N)=ZN^i;
end 
Zss=kron(eye(2), Zs);
spy(Zg)
H=[1 1; 1 -1]/sqrt(2);
Hg=kron(H, kron(eye(N), kron(eye(2), eye(N))));

M=Hg*Zss*Zg;

% A=reshape(M*phi_init, N2, []);
% error= norm(A(1:N, 1:N)-S/2)
% [Ux,r]=qr(x); Ux=r(1,1)*Ux;
% [Uy,r]=qr(y); Uy=r(1, 1)*Uy;
% U1=kron(H,kron(Uy, kron(eye(2), Ux)));
% [q, r]=qr(phi_init);
% U2=q; U2=r(1,1)*U2;
In = eye(N);
K = zeros(4*N^3,N);
Kb = K;
b=rand(N, 1);
b=b/norm(b);
Sb=S*b/(sqrt(2)^(log2(N)+2));

for j = 1:N
    K(:,j) = kron(In(:, j), phi_init);
    Kb(:, j)=kron(b, phi_init);
 %    K2(:,j) = kron(kron(y, x),In(:,j));
end
% for j = 1:N
%     K(:,j) = kron( phi_init, In(:, j));
%     K2(:,j) = kron(kron(y, x),In(:,j));

% end


ZZ=diag(ones(4*N^2-1, 1), 1); ZZ(4*N^2, 1)=1;
Zus=ZZ^(2*N); %effettua lo shift di 2N posizioni in alto

ZZus=eye(4*N^3);
for i=1:N-1
    ZZus(i*4*N^2+1:(i+1)*4*N^2, i*4*N^2+1:(i+1)*4*N^2)=Zus^i;
end
BB=ZZus*kron(eye(N), M)*K;


BBb=ZZus*kron(eye(N), M)*Kb;
Hn=1;
for j=1:log2(N)
    Hn=kron(H, Hn);
end
AB=kron(Hn, eye(4*N^2))*BB;
%AB=kron(eye(4*N^2), Hn)*BB;
ABb=kron(Hn, eye(4*N^2))*BBb;

%Questo realizza un block encoding con alpha= sqrt(2)^(n+2), e servono 3n+2
%qubits aggiuntivi
norm(AB(1:N, 1:N)-S/(2*sqrt(2)^log2(N)))
norm(ABb(1:N, 1)-Sb)
[q, r]=qr(kron(In, phi_init));
U_init=q; U_init=U_init*diag([diag(r); ones(4*N^3-N, 1)]);

tot_wires = 3*n + 2;
circ = qclab.QCircuit(tot_wires);
% e=eye(N);
% eg=eye(4*N^3);
% perm=eye(4*N^3);
% perm(1:N,1:N)=zeros(N,N);
% for i=0:N-1
%     perm=perm+eg(:,i+1)*eg(i*4*N^2+1,:)+eg(:,i*4*N^2+1)*eg(i+1,:);
% end
% perm(1,1)=1;
% spy(perm)
% circ.push_back(qclab.qgates.MatrixGate(0:tot_wires-1,perm))

for i=0:n-1
    circ.push_back(qclab.qgates.SWAP(i,2*n+2+i))
end
Ux=make_circuit(y,n+1);
Ux.asBlock("Uy");
Uy=make_circuit(x,2*n+2);
Uy.asBlock("Uy");

circ.push_back(Ux);

circ.push_back(Uy);
% circ.push_back(qclab.qgates.Hadamard(n-1))
% for i=0:n-1
%     circ.push_back(qclab.qgates.Hadamard(i))
% end
circ.push_back(init_circ(n));


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

% S
% myS
% %%%%%%%%%%%%%%%%%%%%%%%
% uxmat=Ux.matrix
% y

