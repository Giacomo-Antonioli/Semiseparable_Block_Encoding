n=3;
N=2^n;

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
BB=ZZus*kron(eye(N), M);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
circ=BB_gate(n);
mat=circ.matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Look for existing figure
fig = findobj('Type', 'Figure', 'Tag', 'SpyFigure');

if isempty(fig)
    fig = figure('Tag', 'SpyFigure');
else
    figure(fig);   % bring existing figure to foreground
    clf(fig);      % clear previous content
end
subplot(1, 2, 1);
spy(BB);
title('Original');

subplot(1, 2, 2);
spy(mat);
title('Circuit');

norm(abs(BB)-abs(mat))