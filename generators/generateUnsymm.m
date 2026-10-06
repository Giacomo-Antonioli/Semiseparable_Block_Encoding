function [S,u,v,x,y,nu,nv,nx,ny] = generateUnsymm(N)

u=randn(N,1);v=randn(N,1);x=randn(N,1);

y=(u.*v)./x;
nu=norm(u);
nv=norm(v);
nx=norm(x);
ny=norm(y);
% u=u/nu;
% v=v/nv;
% x=x/nx;
% y=y/ny;
S=tril(u*v')+triu(x*y', 1);

end