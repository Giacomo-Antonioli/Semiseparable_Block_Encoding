function [u,v]=uv(s)
% data la matrice s simmetrica, costruisce u e v che definiscono una semiseparabile
% che modella prime riga e colonna
n=size(s,1);
u=s(:,1); v=zeros(1,n);
k=1;
while s(k,1)==0 & k<n k=k+1; end
if k==n & s(n,1)==0 return; else
    v(1)=1; k=n;
    while u(k)==0 & k>1 k=k-1; end
    for j=2:k v(j)=s(k,j)/u(k); end
end
p=sqrt(norm(v)/norm(u)); % bilanciamento di un e vn
u=p*u; v=v/p;