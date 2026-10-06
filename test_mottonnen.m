%TEST_MOTTONEN_STATE_PREP_QCLAB
% Verifies mottonen_state_prep_qclab against random target states:
%   (1) circ.matrix is unitary
%   (2) circ.matrix * e1 matches psi up to global phase
%   (3) cross-checks against the dense-matrix mottonen_state_prep.m
%       implementation (should agree exactly, since both encode the
%       same algorithm)

clear;
clc;

%% ============================================================
% CONFIGURATION
%% ============================================================

n_range    = 1:5;
num_trials = 10;

%% ============================================================
% MAIN VERIFICATION LOOP
%% ============================================================

fprintf('%4s | %10s | %12s | %12s | %14s\n', ...
    'n', 'N', 'unitarity', 'prep err', 'vs dense err');
fprintf('%s\n', repmat('-', 1, 62));

max_unitarity_err = 0;
max_prep_err       = 0;
max_dense_err      = 0;

for n = n_range
    N = 2^n;

    for r = 1:num_trials
        rng(100*n + r);

        psi = randn(N, 1) + 1i * randn(N, 1);
        psi = psi / norm(psi);

        % --- QCLAB circuit ---
        [circuit, global_phase] = mottonen_state_prep_qclab(psi);
        U = get_circuit_matrix(circuit, N);

        % (1) unitarity
        unitarity_err = norm(U' * U - eye(N), 'fro');

        % (2) prepared state matches psi up to global phase
        e1 = zeros(N, 1); e1(1) = 1;
        prepared = U * e1;

        % align phases: find the best-fit global phase between
        % `prepared` and `psi`, then compare
        overlap = psi' * prepared;              % complex number
        phase_fix = overlap / abs(overlap);      % unit-modulus phase alignment
        prep_err = norm(prepared * conj(phase_fix) - psi);

        % also check that the algorithm's own reported global_phase
        % lines up with the phase alignment we just computed
        reported_phase_err = abs(mod(angle(phase_fix) + global_phase + pi, 2*pi) - pi);

        % --- cross-check against the dense-matrix implementation ---
        % Both implementations should encode exactly the same unitary
        % up to a single overall global phase (the dense version applies
        % its own global_phase correction; the QCLAB circuit doesn't).
        circ_dense = mottonen_state_prep_classical(psi);
        M1 = circ_dense.matrix;
        M2 = U;
        k = find(abs(M1(:)) > 1e-8, 1);
        rel_phase = M1(k) / M2(k);
        dense_err = norm(M1 - rel_phase * M2, 'fro');

        max_unitarity_err = max(max_unitarity_err, unitarity_err);
        max_prep_err      = max(max_prep_err, prep_err);
        max_dense_err     = max(max_dense_err, dense_err);

        if unitarity_err > 1e-8
            warning('n=%d trial=%d: unitarity error = %.3e', n, r, unitarity_err);
        end
        if prep_err > 1e-8
            warning('n=%d trial=%d: state prep error = %.3e', n, r, prep_err);
        end
        if reported_phase_err > 1e-6
            warning('n=%d trial=%d: reported global_phase disagrees with actual phase (err=%.3e)', ...
                n, r, reported_phase_err);
        end
        if dense_err > 1e-8
            warning('n=%d trial=%d: mismatch vs dense implementation = %.3e', n, r, dense_err);
        end
    end

    fprintf('%4d | %10d | %12.3e | %12.3e | %14.3e\n', ...
        n, N, max_unitarity_err, max_prep_err, max_dense_err);
end

fprintf('\nAll errors should be at the level of floating-point roundoff (~1e-13 or smaller).\n');

%% ============================================================
% LOCAL FUNCTIONS
%% ============================================================

function U = get_circuit_matrix(circuit, N)
% Returns the N x N unitary matrix implemented by a qclab.QCircuit.
    U = circuit.matrix;
    if ~isequal(size(U), [N, N])
        error('get_circuit_matrix:badSize', ...
            'circuit.matrix returned a %dx%d matrix, expected %dx%d.', ...
            size(U,1), size(U,2), N, N);
    end
end