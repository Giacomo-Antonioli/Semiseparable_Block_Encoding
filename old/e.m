function v = e(k, n)
% E  Return the computational basis vector |k> in dimension 2^n.
%
%   v = e(k, n)
%
%   Inputs:
%       k - index of the basis state (0 ≤ k < 2^n)
%       n - number of qubits
%
%   Output:
%       v - column vector of length 2^n with a 1 at position k+1

    N = 2^n;
    assert(k >= 0 && k < N, 'Index k must satisfy 0 ≤ k < 2^n.');

    v = zeros(N, 1);
    v(k + 1) = 1;   % MATLAB is 1-based indexing
end
