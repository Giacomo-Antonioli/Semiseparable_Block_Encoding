% clear;
% clc;
% N = 16; rng(42);
% x = randn(N,1); x = x/norm(x);
% y = randn(N,1); y = y/norm(y);
% S = tril(x*y') + triu(y*x', 1);
% 
% circ = build_semiseparable_circuit(x, y);
% %circ = build_semiseparable_circuit2(x, y);
% 
% 
% 
% circ.draw
% res=circ.matrix;
% mat=res(1:N,1:N)*2*sqrt(N)
% norm(mat-S)
% 
% 
% %% test_statistics.m
% % Test the circuit statistics methods using a known circuit.
% % This script builds a simple circuit with gates, a barrier, and a subcircuit
% % (block), then checks that getQubitCount, getGateCount, getDepth,
% % and getControlledStats produce the expected values.
% % 
% % clear; clc;
% % 
% % % ------------------------------------------------------------
% % % 1. Build a test circuit
% % % ------------------------------------------------------------
% % 
% % nbQubits = 3;
% % circ = qclab.QCircuit(nbQubits);
% % 
% % % Single‑qubit gates
% % circ.push_back(qclab.qgates.Hadamard(0));
% % circ.push_back(qclab.qgates.PauliX(1));
% % circ.push_back(qclab.qgates.RotationY(2, pi/4));
% % 
% % % Single‑controlled gate: CNOT
% % circ.push_back(qclab.qgates.CNOT(0, 2));
% % 
% % % Multi‑controlled gate: MCX with 2 controls
% % circ.push_back(qclab.qgates.MCX([0, 1], 2));
% % 
% % % Barrier (should be ignored everywhere)
% % circ.barrier();
% % 
% % % Subcircuit (block) with gates
% % block = qclab.QCircuit(nbQubits);
% % block.push_back(qclab.qgates.CRotationZ(1, 2, pi/3));   % single‑controlled rotation
% % block.push_back(qclab.qgates.RotationZ(0, pi/2));
% % block.asBlock('TestBlock');
% % 
% % % Add the block to the main circuit
% % circ.push_back(block);
% % 
% % % Another single‑qubit gate after the block
% % circ.push_back(qclab.qgates.PauliZ(1));
% % 
% % % ------------------------------------------------------------
% % % 2. Manually compute expected values
% % % ------------------------------------------------------------
% % 
% % % Expected gate count (barriers excluded):
% % %   - 3 single‑qubit: Hadamard, PauliX, RotationY
% % %   - 1 CNOT (single‑controlled)
% % %   - 1 MCX (multi‑controlled)
% % %   - 2 gates inside block: CRotationZ (single‑controlled), RotationZ
% % %   - 1 PauliZ
% % % Total = 3 + 1 + 1 + 2 + 1 = 8
% % expected_gates = 8;
% % 
% % % Controlled gates:
% % %   - Single‑controlled: CNOT (1) + CRotationZ (1) = 2
% % %   - Multi‑controlled: MCX (1) = 1
% % expected_single = 2;
% % expected_multi  = 1;
% % expected_total  = expected_single + expected_multi;
% % 
% % % Depth (with unit‑time gates, barriers ignored):
% % %   Layer 1: H(0), X(1), RY(2)   → parallel (different qubits)
% % %   Layer 2: CNOT(0,2)
% % %   Layer 3: MCX(0,1,2)
% % %   Barrier ignored
% % %   Layer 4: CRotationZ(1,2) and RotationZ(0) can be parallel
% % %            because they act on disjoint qubits (1,2) and (0)
% % %   Layer 5: PauliZ(1) (after block)
% % expected_depth = 5;
% % 
% % % Number of qubits
% % expected_qubits = nbQubits;
% % 
% % % ------------------------------------------------------------
% % % 3. Run statistics methods
% % % ------------------------------------------------------------
% % circ.draw
% % fprintf('=== Circuit Statistics Test ===\n\n');
% % 
% % % getQubitCount
% % qubits = circ.getQubitCount();
% % fprintf('Qubit count: %d (expected %d)  %s\n', qubits, expected_qubits, ...
% %     check(qubits, expected_qubits));
% % 
% % % getGateCount
% % gates = circ.getGateCount();
% % fprintf('Gate count : %d (expected %d)  %s\n', gates, expected_gates, ...
% %     check(gates, expected_gates));
% % 
% % % getDepth
% % depth = circ.getDepth();
% % fprintf('Depth      : %d (expected %d)  %s\n', depth, expected_depth, ...
% %     check(depth, expected_depth));
% % 
% % % getControlledStats
% % ctrl = circ.getControlledStats();
% % fprintf('Controlled gates:\n');
% % fprintf('  Single   : %d (expected %d)  %s\n', ctrl.single, expected_single, ...
% %     check(ctrl.single, expected_single));
% % fprintf('  Multi    : %d (expected %d)  %s\n', ctrl.multi, expected_multi, ...
% %     check(ctrl.multi, expected_multi));
% % fprintf('  Total    : %d (expected %d)  %s\n', ctrl.total, expected_total, ...
% %     check(ctrl.total, expected_total));
% % 
% % % ------------------------------------------------------------
% % % 4. Optional: inspect flattened gate list
% % % ------------------------------------------------------------
% % fprintf('\n=== Flattened gate list (objectsFlattened) ===\n');
% % flat = circ.objectsFlattened();
% % for i = 1:length(flat)
% %     gate = flat(i);
% %     % Class name without package prefix
% %     cls = class(gate);
% %     if startsWith(cls, 'qclab.qgates.')
% %         name = cls(13:end);
% %     else
% %         name = cls;
% %     end
% %     % Qubit info
% %     if isprop(gate, 'qubits')
% %         qstr = sprintf('qubits=[%s]', num2str(gate.qubits));
% %     elseif isprop(gate, 'target') && isprop(gate, 'controls')
% %         qstr = sprintf('controls=[%s], target=%d', num2str(gate.controls), gate.target);
% %     elseif isprop(gate, 'qubit')
% %         qstr = sprintf('qubit=%d', gate.qubit);
% %     else
% %         qstr = '';
% %     end
% %     fprintf('Gate %2d: %-20s %s\n', i, name, qstr);
% % end
% % fprintf('Total gates in flattened list: %d\n', length(flat));
% % 
% % % Verify that barriers are present in the flattened list (they are),
% % % but excluded from getGateCount.
% % no_barrier_count = 0;
% % for i = 1:length(flat)
% %     if ~isa(flat(i), 'qclab.qgates.Barrier')
% %         no_barrier_count = no_barrier_count + 1;
% %     end
% % end
% % fprintf('Gates excluding barriers: %d (matches getGateCount: %d)\n', ...
% %     no_barrier_count, circ.getGateCount());
% % 
% % % ------------------------------------------------------------
% % % Helper function for PASS/FAIL
% % % ------------------------------------------------------------
% % function [result] = check(actual, expected)
% %     if actual == expected
% %         result = 'PASS';
% %     else
% %         result = sprintf('FAIL (diff=%d)', actual - expected);
% %     end
% % end

clear; clc;

%% 1. Build the circuit
N = 4; rng(42);
x = randn(N,1); x = x/norm(x);
y = randn(N,1); y = y/norm(y);
S = tril(x*y') + triu(y*x', 1);

%circ = build_semiseparable_circuit(x, y);
circ = build_semiseparable_circuit2(x, y);

%% 2. Sanity check: does the circuit reproduce S?
circ.draw
res = circ.matrix;
mat = res(1:N,1:N)*2*sqrt(N);
fprintf('Reconstruction error ||mat - S|| = %.3e\n\n', norm(mat-S));

%% 3. Circuit statistics
fprintf('=== Circuit Statistics ===\n\n');

qubits = circ.getQubitCount();
fprintf('Qubit count       : %d\n', qubits);

gates = circ.getGateCount();
fprintf('Gate count        : %d\n', gates);

depth = circ.getDepth();
fprintf('Depth             : %d\n', depth);

ctrl = circ.getControlledStats();
fprintf('Controlled gates  :\n');
fprintf('  Single-controlled: %d\n', ctrl.single);
fprintf('  Multi-controlled : %d\n', ctrl.multi);
fprintf('  Total controlled : %d\n', ctrl.total);

%% 4. Optional: inspect flattened gate list
fprintf('\n=== Flattened gate list (objectsFlattened) ===\n');
flat = circ.objectsFlattened();
for i = 1:length(flat)
    gate = flat(i);
    cls = class(gate);
    if startsWith(cls, 'qclab.qgates.')
        name = cls(13:end);
    else
        name = cls;
    end

    if ismethod(gate, 'qubits')
        qstr = sprintf('qubits=[%s]', num2str(gate.qubits));
    elseif ismethod(gate, 'target') && ismethod(gate, 'controls')
        qstr = sprintf('controls=[%s], target=%d', num2str(gate.controls), gate.target);
    elseif ismethod(gate, 'target') && ismethod(gate, 'control')
        qstr = sprintf('control=%d, target=%d', gate.control, gate.target);
    elseif ismethod(gate, 'qubit')
        qstr = sprintf('qubit=%d', gate.qubit);
    else
        qstr = '';
    end
    fprintf('Gate %3d: %-20s %s\n', i, name, qstr);
end
fprintf('Total gates in flattened list: %d\n', length(flat));

no_barrier_count = 0;
for i = 1:length(flat)
    if ~(isa(flat(i), 'qclab.Barrier') || isa(flat(i), 'qclab.qgates.Barrier'))
        no_barrier_count = no_barrier_count + 1;
    end
end
fprintf('Gates excluding barriers: %d (matches getGateCount: %d)\n', ...
    no_barrier_count, circ.getGateCount());