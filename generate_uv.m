function [u,v] = generate_uv(n, type)

    N = 2^n;

    switch type

        %% =========================
        % TYPE 1: RANDOM
        %% =========================
        case 1
            u = randn(N,1); u = u/norm(u);
            v = randn(N,1); v = v/norm(v);

        %% =========================
        % TYPE 2: TRIDIAGONAL + uv(s)
        %% =========================
        case 2
            T = diag(randn(N,1)) + ...
                diag(randn(N-1,1),1) + ...
                diag(randn(N-1,1),-1);

            T = T + 1e-2*eye(N);

            x = randn(N,1);
            y = randn(N,1);

            s = (T\x) * (T\y)';   % structured rank-1
            s = (s + s')/2;

            [u,v] = uv(s);
            v=v';
        %% =========================
        % TYPE 3: EXPONENTIAL OU
        %% =========================
        case 3
            gamma = cumsum(randn(N,1));

            u = exp(gamma);
            v = exp(-gamma);

            u = u/norm(u);
            v = v/norm(v);

        otherwise
            error('type must be 1, 2, or 3');
    end
end