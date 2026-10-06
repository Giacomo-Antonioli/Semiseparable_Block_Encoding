function [u, v] = gen_ornstein_uhlenbeck(n, theta)
    % GEN_ORNSTEIN_UHLENBECK Generators of the covariance matrix of a 
    % stationary Ornstein-Uhlenbeck (Gauss-Markov) process sampled on a 
    % uniform grid of [0,1].
    %
    % C_ij = exp(-theta |t_i - t_j|) = u_i v_j for i<=j, with
    %     u_i = exp(theta * t_i),  v_j = exp(-theta * t_j).
    % Returned unit-normalized.
    
    % Set default theta if not provided
    if nargin < 2
        theta = 3.0;
    end

    N = 2^n;
    t = linspace(0.0, 1.0, N); % Row vector by default
    
    u = exp(theta * t);
    v = exp(-theta * t);
    
    % Unit-normalize using the Euclidean norm
    u = u / norm(u);
    v = v / norm(v);
end