n=3;
N=2^n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
circ=Zbi_gate(n);
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
spy(Zbi);
title('Original');

subplot(1, 2, 2);
spy(mat);
title('Circuit');

norm(abs(Zbi)-abs(mat))