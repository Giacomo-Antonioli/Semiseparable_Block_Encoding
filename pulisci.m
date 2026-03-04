function[a]=pulisci(a, tol)
[n,m]=size(a);
for i=1:n
    for k=1:m
        if(abs(a(i,k))<tol)
            a(i,k)=0;
        end
    end
end