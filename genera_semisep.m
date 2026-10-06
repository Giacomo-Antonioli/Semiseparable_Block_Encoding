function [S, u, v, p, q]=genera_semisep(n)
u= randn(n, 1); v=randn(n, 1); u(n-1:n)=0.5*randn(2, 1); v(n-1:n)=0.5*randn(2, 1);
u= u/norm(u); v=v/norm(v);
d= u.*v;
t= zeros(n, 1);
w=abs(d);
l=ones(n, 1); % rappresenta exp(2t_i), le prime n-2 componenti sono 1
s= 1-sum(w(3:n));
delta= (s^2-w(1)^2+w(2)^2)^2-4*s^2*w(2)^2;
if delta<0
    error('generatori non buoni');
end
l(2)=(s^2-w(1)^2+w(2)^2+sqrt(delta))/(2*w(2)*s);
l(1)= (s-w(2)*l(2))/w(1);

%t=log(l)/2;
p= sqrt(abs(d)).*sqrt(l);
q= sign(d).*sqrt(abs(d)).*(1./sqrt(l));
S=tril(u*v')+triu(p*q', 1);
if abs(norm(p)-1)>10^-8 || abs(norm(q)-1)>10^-8
    error('problema di stabilità')
end