%% Parameters
n = 2;              % choose n
N = 2^n;            % size of first identity
M = 2^(n+1);        % size of shift operator

%% Define Z_{2^(n+1)} as cyclic shift
Z = zeros(M);
for k = 1:M
    Z(k, mod(k,M)+1) = 1;
end

%% Construct sum_j |j><j| ⊗ Z^{-j}
sum_op = zeros(M*M, M*M);  % preallocate
for j = 0:M-1
    projector = zeros(M);
    projector(j+1,j+1) = 1;   % |j><j|
    sum_op = sum_op + kron(projector, Z^(-j));
end

%% Final Zg operator
Zg = kron(eye(N), sum_op);

%% Optional: display size
disp(['Size of Zg: ', num2str(size(Zg,1)), ' x ', num2str(size(Zg,2))]);