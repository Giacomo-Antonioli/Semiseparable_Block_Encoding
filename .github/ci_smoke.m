function ci_smoke()
%CI_SMOKE Build the symmetric block encoding and check it against S.
    setup_paths();
    assert(~isempty(which('qclab.qgates.CSWAP')), 'bundled QCLAB fork not on the path');

    for N = [2 4 8]
        rng(42);
        u = randn(N,1);  u = u / norm(u);
        v = randn(N,1);  v = v / norm(v);
        S = tril(u*v') + triu(v*u', 1);

        circ = build_semiseparable_circuit(u, v);
        U = circ.matrix;

        unitarity = norm(U' * U - eye(size(U)), 2);
        block     = norm(2*sqrt(N)*U(1:N, 1:N) - S, 2);
        fprintf('N=%d  unitarity %.2e  block %.2e  gates %d  depth %d\n', ...
                N, unitarity, block, circ.getGateCount, circ.getDepth);
        assert(unitarity < 1e-10 && block < 1e-10, 'block encoding check failed for N=%d', N);
    end
end
