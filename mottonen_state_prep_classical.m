function circ = mottonen_state_prep_classical(psi)
%MOTTONEN_STATE_PREP  Möttönen et al. state preparation via uniformly
%controlled (multiplexed) rotations.
%
%   circ = MOTTONEN_STATE_PREP(psi) builds the explicit unitary
%   U (N x N, N = 2^n) such that U*e1 = psi, where e1 = [1;0;...;0]
%   is the all-zero computational basis state |0...0>, using the
%   algorithm of:
%
%       Möttönen, Vartiainen, Bergholm, Salomaa,
%       "Transformation of quantum states using uniformly controlled
%       rotations", Quantum Inf. Comput. 5(6):467-473 (2005),
%       arXiv:quant-ph/0407010.
%
%   INPUT
%       psi : complex column vector of length N = 2^n. Need not be
%             exactly normalized (it is renormalized internally; a
%             warning is issued if it deviates from unit norm by more
%             than 1e-8).
%
%   OUTPUT (struct circ)
%       circ.matrix      : N x N unitary such that circ.matrix(:,1) == psi
%       circ.n           : number of qubits
%       circ.ry_angles   : 1 x n cell array, ry_angles{k} is the vector
%                          of 2^(k-1) multiplexed-Ry rotation angles
%                          applied to qubit k
%       circ.rz_angles   : 1 x n cell array, same structure for Rz
%       circ.global_phase: residual global phase applied at the end
%
%   QUBIT / BASIS CONVENTION
%       Qubit 1 is the most significant bit of the basis-state index.
%       Gates are applied in the order: Ry on qubit 1 (uncontrolled),
%       Ry on qubit 2 (controlled by qubit 1), ..., Ry on qubit n
%       (controlled by qubits 1..n-1), then the same pattern for Rz,
%       then a global phase correction so that circ.matrix(:,1) equals
%       psi exactly (not just up to global phase).
%
%   NOTE
%       This builds the *exact* unitary directly via block-diagonal
%       construction (kron/blkdiag) for verification/simulation
%       purposes. It does not perform the further decomposition of
%       each multiplexed rotation into a Gray-code sequence of
%       single-qubit rotations and CNOTs that a real gate-count-
%       efficient circuit would use.

    psi = psi(:);
    N = numel(psi);
    n = round(log2(N));

    if 2^n ~= N
        error('mottonen_state_prep:badLength', ...
            'Length of psi must be a power of two (got %d).', N);
    end

    nrm = norm(psi);
    if abs(nrm - 1) > 1e-8
        warning('mottonen_state_prep:notNormalized', ...
            '||psi|| = %.6g != 1; renormalizing.', nrm);
    end
    psi = psi / nrm;

    %% ------------------------------------------------------------
    % Build the magnitude/phase tree, level by level, from the
    % leaves (level n, the amplitudes themselves) up to the root
    % (level 0, a single node = the overall norm / overall phase).
    %% ------------------------------------------------------------

    r_k   = abs(psi);      % level-n magnitudes, length N
    phi_k = angle(psi);    % level-n phases, length N

    ry_angles = cell(1, n);
    rz_angles = cell(1, n);

    for k = n:-1:1
        R = reshape(r_k, 2, []);       % row1 = even (bit=0), row2 = odd (bit=1)
        Phi = reshape(phi_k, 2, []);

        ry_angles{k} = (2 * atan2(R(2, :), R(1, :))).';   % length 2^(k-1)
        rz_angles{k} = (Phi(2, :) - Phi(1, :)).';         % length 2^(k-1)

        r_k   = hypot(R(1, :), R(2, :)).';                % length 2^(k-1)
        phi_k = ((Phi(1, :) + Phi(2, :)) / 2).';          % length 2^(k-1)
    end

    global_phase = phi_k(1);   % what's left at the root: r_k(1) should be ~1

    %% ------------------------------------------------------------
    % Assemble the full unitary from the per-level multiplexed gates.
    %% ------------------------------------------------------------

    U = eye(N);
    for k = 1:n
        U = build_multiplexed_gate(n, k, ry_angles{k}, 'ry') * U;
    end
    for k = 1:n
        U = build_multiplexed_gate(n, k, rz_angles{k}, 'rz') * U;
    end
    U = exp(1i * global_phase) * U;

    circ.matrix       = U;
    circ.n            = n;
    circ.ry_angles    = ry_angles;
    circ.rz_angles    = rz_angles;
    circ.global_phase = global_phase;
end

%% ==================================================================
% LOCAL FUNCTIONS
%% ==================================================================

function M = build_multiplexed_gate(n, k, angles, gate_type)
% Builds the N x N (N = 2^n) unitary for a uniformly controlled
% rotation targeting qubit k, controlled by qubits 1..k-1, acting as
% identity on qubits k+1..n. angles must have length 2^(k-1); angles(c+1)
% is the rotation angle used when qubits 1..k-1 are in basis state c
% (0-indexed).

    num_blocks     = numel(angles);         % = 2^(k-1)
    spectator_dim  = 2^(n - k);
    I_spectator    = eye(spectator_dim);

    blocks = cell(1, num_blocks);
    for c = 1:num_blocks
        theta = angles(c);
        switch gate_type
            case 'ry'
                g = [cos(theta/2), -sin(theta/2); ...
                     sin(theta/2),  cos(theta/2)];
            case 'rz'
                g = diag([exp(-1i*theta/2), exp(1i*theta/2)]);
            otherwise
                error('build_multiplexed_gate:badType', ...
                    'Unknown gate_type "%s"', gate_type);
        end
        blocks{c} = kron(g, I_spectator);
    end

    M = blkdiag(blocks{:});
end