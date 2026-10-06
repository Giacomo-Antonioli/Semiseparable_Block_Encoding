function [u, v] = local_generate_uv(n, type)

    N = 2^n;

    switch type

        case 1
            u = randn(N, 1); u = u / norm(u);
            v = randn(N, 1); v = v / norm(v);

        case 2
            % T = diag(randn(N, 1)) + ...
            %     diag(randn(N-1, 1), 1) + ...
            %     diag(randn(N-1, 1), -1);
            % T=T+T'+2*eye(N);
            % S = inv(T);
            % [u, v] = uv(S);
            % v = v';
            % u = u / norm(u);
            % v = v / norm(v);
             [u, v, ~, ~] =ss_generators_from_tridiag(N, 123, 2);

        case 3
            [u, v] = gen_ornstein_uhlenbeck(n, 3.0);
            u=u';
            v=v';
        otherwise
            error('local_generate_uv:badType', 'Unknown type %d', type);
    end
end