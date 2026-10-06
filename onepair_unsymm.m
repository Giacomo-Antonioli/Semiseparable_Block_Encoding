function [S,u,v,x,y] = onepair_unsymm(N)
% One-pair (generator-representable) NON-symmetric semiseparable matrix.
% S_ij = x_i y_j  (i<=j),  u_j v_i (i>=j), with diagonal compatibility x_i y_i = u_i v_i.
    
    u = randn(N,1);  v = randn(N,1);  x = randn(N,1);
    x(x==0) = 1;                    % avoid division by zero
    y = (u.*v)./x;                  % enforce  x_i y_i = u_i v_i  on the diagonal
    % unit-norm generators (optional, to match the theorem hypotheses)
    u=u/norm(u); v=v/norm(v); x=x/norm(x); y=y/norm(y);
    % rebuild (renormalization breaks compatibility, so re-enforce y afterwards)
    y = (u.*v)./x;                  % re-enforce after normalizing u,v,x
    y=y/norm(y);
    [I,J] = ndgrid(1:N,1:N);
    S = zeros(N);
    upper = I<=J;  lower = ~upper;
    S(upper) = x(I(upper)).*y(J(upper));
    S(lower) = u(J(lower)).*v(I(lower));
end