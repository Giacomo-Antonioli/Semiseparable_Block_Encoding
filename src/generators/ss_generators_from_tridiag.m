function [u, v, scal, A] = ss_generators_from_tridiag(N, seed, ddom)
%SS_GENERATORS_FROM_TRIDIAG  Generators of the semiseparable inverse of a
%   random irreducible symmetric tridiagonal matrix.
%
%   [u,v,scal,A] = ss_generators_from_tridiag(N,seed,ddom)
%
%   Builds a random symmetric tridiagonal A (irreducible: all off-diagonal
%   entries bounded away from zero; diagonally dominant with strength DDOM to
%   guarantee invertibility and control the conditioning), and returns the
%   generators u,v of its inverse, which is a one-pair semiseparable matrix:
%
%        inv(A)(i,j) = u(min(i,j)) * v(max(i,j)) / scal
%
%   u,v are returned UNIT-NORM (as required by the hypotheses of the block-
%   encoding theorem); SCAL is the scalar linking S(u,v) to inv(A):
%
%        S(u,v) = scal * inv(A),        inv(A) = S(u,v)/scal.
%
%   STABILITY.  The classical formula for the inverse of a tridiagonal matrix
%   involves products of leading/trailing principal minors, which under- and
%   overflow exponentially in N and must NOT be used.  Instead we exploit the
%   fact that the first and last columns of inv(A) are already proportional to
%   the two generators:
%        A \ e_1   = u(1) * v          (first column of inv(A))
%        A \ e_N   = v(N) * u          (last column of inv(A))
%   Each is obtained by one backward-stable tridiagonal solve in O(N) flops.
%   The coupling constant c = u(1)v(N) = inv(A)(1,N) is computed once and used
%   only as a SCALAR (never to rescale a vector), so the exponentially small
%   value of c never contaminates the generators.

if nargin < 2  isempty(seed), seed = 0;    end
if nargin < 3  isempty(ddom), ddom = 1.0;  end   % diagonal-dominance margin

rng(seed);

% ---- random irreducible symmetric tridiagonal ---------------------------
b = randn(N-1,1);
b(abs(b) < 0.1) = 0.1;                 % irreducibility: b_i ~= 0, well separated
a = randn(N,1) + abs([b;0]) + abs([0;b]) + ddom;   % diagonal dominance
A = spdiags([[b;0] a [0;b]], -1:1, N, N);          % sparse tridiagonal

% ---- stable extraction of the generators --------------------------------
e1 = sparse(1,1,1,N,1);
eN = sparse(N,1,1,N,1);

p = A \ e1;        % = u(1) * v      (first column of inv(A))
q = A \ eN;        % = v(N) * u      (last  column of inv(A))
p = full(p);  q = full(q);

c = p(N);          % = u(1)*v(N) = inv(A)(1,N)   (equals q(1); see check below)

% consistency check (should agree to roundoff)
if abs(p(N) - q(1)) > 1e-8 * max(abs(p(N)), eps)
    warning('ss_generators:asym', ...
        'p(N) and q(1) differ (%g vs %g): A may be ill-conditioned.', p(N), q(1));
end

nu = norm(q);  nv = norm(p);
u  = q / nu;       % unit-norm generator (u)
v  = p / nv;       % unit-norm generator (v)

scal = c / (nu*nv);   % S(u,v) = scal * inv(A)
end