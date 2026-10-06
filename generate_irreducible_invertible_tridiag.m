function A = generate_irreducible_invertible_tridiag(N)
    % N: size of the square matrix
    
    % 1. Generate off-diagonals strictly away from 0 to ensure irreducibility.
    % We use random values in [0.5, 1.5] to guarantee they are never zero.
    l = 0.5 + rand(N-1, 1); % Subdiagonal
    u = 0.5 + rand(N-1, 1); % Superdiagonal
    
    % 2. Compute the main diagonal to ensure strict diagonal dominance
    d = zeros(N, 1);
    for i = 1:N
        if i == 1
            sum_off = abs(u(1));
        elseif i == N
            sum_off = abs(l(N-1));
        else
            sum_off = abs(l(i-1)) + abs(u(i));
        end
        
        % Add a positive buffer (e.g., +1.0) to make it strictly dominant
        d(i) = sum_off + 1.0; 
    end
    
    % 3. Assemble the tridiagonal matrix
    A = diag(d) + diag(l, -1) + diag(u, 1);
end