function circuit = mottonen_state_prep_qclab(psi, offset)
%MOTTONEN_STATE_PREP_QCLAB  Möttönen et al. state preparation as a QCLAB
%circuit (Ry/Rz uniformly-controlled rotations + CNOTs + a global-phase
%correction).
%
%   circuit = MOTTONEN_STATE_PREP_QCLAB(psi) returns a qclab.QCircuit on
%   n = log2(length(psi)) qubits such that applying the circuit to
%   |0...0> reproduces psi exactly:
%
%       circuit.matrix * e1  ==  psi / norm(psi)
%
%   where e1 = [1;0;...;0]. The global phase needed to match psi exactly
%   (not just up to phase) is applied as an explicit gate sequence on
%   qubit 0 -- see "GLOBAL PHASE" below.
%
%   circuit = MOTTONEN_STATE_PREP_QCLAB(psi, offset) uses an offset so
%   that the prepared state is encoded on qubits offset:offset+n-1,
%   with the rest of the (n+offset)-qubit register untouched. All gate
%   targets and controls are shifted by offset.
%
%   QUBIT CONVENTION (QCLAB): qubit 0 is the most significant qubit
%   (drawn on top). This matches the "qubit 1 = MSB" convention used to
%   derive the rotation angles, shifted by one (math qubit k -> QCLAB
%   qubit k-1).
%
%   Reference: Möttönen, Vartiainen, Bergholm, Salomaa, "Transformation
%   of quantum states using uniformly controlled rotations",
%   arXiv:quant-ph/0407010.
%
%   Each uniformly-controlled (multiplexed) rotation with target qubit t
%   and s controls is built recursively: split the 2^s angles into the
%   top-control-qubit's two halves, take the elementwise average
%   (a_plus) and half-difference (a_minus), recurse on the remaining
%   s-1 controls for each half, and sandwich as:
%
%       [recurse on a_minus]  CNOT  [recurse on a_plus]  CNOT
%
%   (applied left-to-right in time: a_plus block first). This uses
%   2^(s+1) - 2 CNOTs per multiplexed rotation -- correct, though not
%   the CNOT-count-optimal Gray-code version.
%
%   GLOBAL PHASE
%   RY, RZ, and CNOT all have determinant +-1 when embedded in the full
%   N x N space, so no circuit built purely from them can ever realize
%   a generic global phase e^{i*phi}*I_N (determinant e^{i*phi*N}). The
%   Phase gate P(phi) = diag(1, e^{i*phi}) breaks this: sandwiching it
%   with PauliX on the same qubit,
%
%       P(phi), X, P(phi), X
%
%   gives diag(1,e^{i*phi}) * diag(e^{i*phi},1) = e^{i*phi} * I_2 on
%   that qubit, i.e. e^{i*phi} * I_N on the whole register once
%   tensored with identity on every other (untouched) qubit. This is
%   appended on qubit `offset` at the end of the circuit, using the phi0
%   computed by the magnitude/phase tree recursion below.

    if nargin < 2
        offset = 0;
    end
    % Validate offset
    if ~isscalar(offset) || offset < 0 || offset ~= floor(offset)
        error('mottonen_state_prep_qclab:badOffset', ...
            'offset must be a non‑negative integer.');
    end

    psi = psi(:);
    N = numel(psi);
    n = round(log2(N));

    if 2^n ~= N
        error('mottonen_state_prep_qclab:badLength', ...
            'Length of psi must be a power of two (got %d).', N);
    end

    nrm = norm(psi);
    if abs(nrm - 1) > 1e-8
        warning('mottonen_state_prep_qclab:notNormalized', ...
            '||psi|| = %.6g != 1; renormalizing.', nrm);
    end
    psi = psi / nrm;

    %% ------------------------------------------------------------
    % Magnitude/phase tree (identical to the dense-matrix version):
    % level k (k = n downto 1) yields the 2^(k-1) angles for the
    % multiplexed rotation that targets qubit k (math, 1-indexed).
    %% ------------------------------------------------------------

    r_k   = abs(psi);
    phi_k = angle(psi);

    ry_angles = cell(1, n);
    rz_angles = cell(1, n);

    for k = n:-1:1
        R   = reshape(r_k,   2, []);
        Phi = reshape(phi_k, 2, []);

        ry_angles{k} = (2 * atan2(R(2, :), R(1, :))).';
        rz_angles{k} = (Phi(2, :) - Phi(1, :)).';

        r_k   = hypot(R(1, :), R(2, :)).';
        phi_k = ((Phi(1, :) + Phi(2, :)) / 2).';
    end

    global_phase = phi_k(1);

    %% ------------------------------------------------------------
    % Build the circuit. Math qubit k -> QCLAB qubit (k-1)+offset.
    % Controls for level-k are math qubits 1..k-1, i.e.
    % QCLAB qubits offset .. offset+k-2, in that order.
    %% ------------------------------------------------------------

    circuit = qclab.QCircuit(n + offset);

    for k = 1:n
        target   = (k - 1) + offset;
        controls = (0:(k - 2)) + offset;   % empty when k == 1
        add_multiplexed_rotation(circuit, 'ry', target, controls, ry_angles{k});
    end
    for k = 1:n
        target   = (k - 1) + offset;
        controls = (0:(k - 2)) + offset;
        add_multiplexed_rotation(circuit, 'rz', target, controls, rz_angles{k});
    end

    % --- mandatory global-phase correction (see header) ---
    circuit.push_back(qclab.qgates.Phase(offset, global_phase));
    circuit.push_back(qclab.qgates.PauliX(offset));
    circuit.push_back(qclab.qgates.Phase(offset, global_phase));
    circuit.push_back(qclab.qgates.PauliX(offset));
end

%% ==================================================================
% LOCAL FUNCTIONS
%% ==================================================================

function add_multiplexed_rotation(circuit, gate_type, target, controls, angles)
% Appends the uniformly-controlled rotation (target, controls, angles)
% to circuit, recursively, using the sum/difference + CNOT-sandwich
% decomposition. controls(1) is the most-significant control qubit,
% matching angles' bit-pattern indexing (controls(1)'s bit is the MSB
% of the 0-indexed position into angles). All qubit indices are already
% shifted by offset.

    s = numel(controls);

    if s == 0
        push_rotation(circuit, gate_type, target, angles(1));
        return;
    end

    half = numel(angles) / 2;
    a = angles(1:half);       % top control (controls(1)) == 0
    b = angles(half+1:end);   % top control (controls(1)) == 1

    a_plus  = (a + b) / 2;
    a_minus = (a - b) / 2;

    c1   = controls(1);
    rest = controls(2:end);

    add_multiplexed_rotation(circuit, gate_type, target, rest, a_plus);
    circuit.push_back(qclab.qgates.CNOT(c1, target));
    add_multiplexed_rotation(circuit, gate_type, target, rest, a_minus);
    circuit.push_back(qclab.qgates.CNOT(c1, target));
end

function push_rotation(circuit, gate_type, target, theta)
    switch gate_type
        case 'ry'
            circuit.push_back(qclab.qgates.RotationY(target, theta));
        case 'rz'
            circuit.push_back(qclab.qgates.RotationZ(target, theta));
        otherwise
            error('push_rotation:badType', 'Unknown gate_type "%s"', gate_type);
    end
end